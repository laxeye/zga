#!/usr/bin/env python3
import argparse
import os.path
import logging
import sys
import shutil
import subprocess
import re
import hashlib
from pathlib import Path
from Bio import SeqIO
from zga import __version__
from zga.assemblers import assemble
from zga.essential import run_external, create_subdir
from zga.preprocessing import read_processing


def parse_args():
	'''Returns argparse.Namespace'''
	parser = argparse.ArgumentParser(prog="zga",
		description=f"ZGA genome assembly and annotation pipeline ver. {__version__}")
	zga_steps = [
		"readqc", "processing", "assembling",
		"polishing", "check_genome", "annotation"]
	general_args = parser.add_argument_group(
		title="General options", description="")
	general_args.add_argument("-s", "--first-step",
		choices=zga_steps, default="readqc",
		help="First step of the pipeline. Default: readqc")
	general_args.add_argument("-l", "--last-step",
		choices=zga_steps, default="annotation",
		help="Last step of the pipeline. Default: annotation")
	general_args.add_argument("-o", "--output-dir", required=True,
		help="Output directory")
	general_args.add_argument("--force", action="store_true",
		help="Overwrite output directory if exists")
	general_args.add_argument("-t", "--threads", type=int, default=4,
		help="Number of CPU threads to use (where possible), default: 4.")
	general_args.add_argument("-m", "--memory-limit", type=int, default=8,
		help="Memory limit in GB, default: 8.")
	general_args.add_argument("--genus", default="Unknown",
		help="Genus name if known, default: 'Unknown'.")
	general_args.add_argument("--species", default="sp.",
		help="Species name if known, , default: 'sp.'")
	general_args.add_argument("--strain",
		help="Strain name if known")
	general_args.add_argument("--transparent", action="store_true",
		help="Show output from tools inside the pipeline")
	general_args.add_argument("--domain",
		default="Bacteria", choices=['Archaea', 'Bacteria'],
		help="Prokaryotic domain: Bacteria or Archaea, default: Bacteria.")
	parser.add_argument("-V", "--version",
		action="version", version=f"ZGA ver. {__version__}")

	# Input
	input_args = parser.add_argument_group(title="Input files and options",
		description="Sequencing reads should be in FASTQ format and may be GZipped. "
		+ "Multiple libraries should be provided as space-separated list. "
		+ "If some type of short reads are partialyy available use n/a. "
		+ "e.g. -1 Monday.R1.fq Friday.R1.fq -2 Monday.R2.fq Friday.R2.fq "
		+ "--pe-merged n/a Friday.merged.fq")
	input_args.add_argument("-1", "--pe-1", nargs='+',
		help="FASTQ file(s) with first (left) paired-end reads. "
		+ "Space-separated if multiple.")
	input_args.add_argument("-2", "--pe-2", nargs='+',
		help="FASTQ file(s) with second (right) paired-end reads. "
		+ "Space-separated if multiple.")
	input_args.add_argument("--pe-merged", nargs='+',
		help="FASTQ file(s) with merged overlapped paired-end reads")
	input_args.add_argument("-S", "--single-end", nargs='+',
		help="FASTQ file(s) with unpaired or single-end reads")
	input_args.add_argument("--mp-1", nargs='+',
		help="Mate pair forward reads. SPAdes only")
	input_args.add_argument("--mp-2", nargs='+',
		help="Mate pair forward reads. SPAdes only")
	input_args.add_argument("--pacbio", nargs='+',
		help="PacBio reads. Space-separated if multiple.")
	input_args.add_argument("--nanopore", nargs='+',
		help="Nanopore reads. Space-separated if multiple.")

	# Read processing
	reads_args = parser.add_argument_group(title="Read processing settings")
	reads_args.add_argument("-q", "--quality-cutoff", type=int, default=18,
		help="Base quality cutoff for short reads, default: 18")
	reads_args.add_argument("--adapters",
		help="Adapter sequences for short reads trimming (FASTA). "
		+ "By default Illumina and BGI adapter sequences are used.")
	reads_args.add_argument("--filter-by-tile", action="store_true",
		help="Filter short reads (Illumina only!) "
		+ "based on positional quality over a flowcell.")
	reads_args.add_argument("--min-short-read-length", type=int, default=33,
		help="Minimum short read length to keep after quality trimming.")
	reads_args.add_argument("--entropy-cutoff", type=float, default=-1,
		help="Set between 0 and 1 to filter reads with entropy below "
		+ "that value. Higher is more stringent. Default = -1, filtering disabled.")
	reads_args.add_argument("--bbduk-k", type=int, default=19,
		help="Kmer length used for finding contaminants with BBduk.")
	reads_args.add_argument("--bbduk-extra", nargs='*',
		help="Extra options for BBduk. Should be space-separated.")
	reads_args.add_argument("--tadpole-correct", action="store_true",
		help="Perform error correction of short reads with tadpole.sh from BBtools."
		+ "SPAdes correction may be disabled with \"--no-spades-correction\".")
	reads_args.add_argument("--bbmerge", action="store_true",
		help="Merge overlapped paired-end reads with BBMerge.")
	reads_args.add_argument("--bbmerge-extra", nargs='*',
		help="Extra options for BBMerge. Should be space-separated.")
	reads_args.add_argument("--normalize-kmer-cov", type=int,
		help="Normalize read depth based on kmer counts to arbitrary value.")
	reads_args.add_argument("--calculate-genome-size", action="store_true",
		help="Estimate genome size with mash.")
	reads_args.add_argument("--genome-size-estimation", type=int,
		help="Genome size in bp (no K/M suffix supported) for Flye assembler, if known.")
	reads_args.add_argument("--mash-kmer-copies", type=int, default=10,
		help="Minimum copies of each k-mer to include in size estimation, default: 10.")
	# Mate pair read processing
	reads_args.add_argument("--use-unknown-mp", action="store_true",
		help="NxTrim: Include reads that are probably mate pairs "
		+ "(default: only known mate pairs used)")
	reads_args.add_argument("--no-nxtrim", action="store_true",
		help="Don't process mate-pair reads with NxTrim. Usefull for preprocessed reads")

	asly_args = parser.add_argument_group(title="Assembly settings")
	asly_args.add_argument("-a", "--assembler",
		default="unicycler", choices=["spades", "unicycler", "flye", "megahit"],
		help="Assembler: Unicycler (default; better quality), "
		+ "SPAdes (faster, may use mate-pair reads), "
		+ "Flye (long reads only with short-reads polishing), "
		+ "MEGAHIT (short reads only).")
	asly_args.add_argument("--no-spades-correction", action="store_true",
		help="Disable short read correction by SPAdes.")
	# Spades options
	asly_args.add_argument("--use-scaffolds", action="store_true",
		help="SPAdes: Use assembled scaffolds. Contigs are used by default.")
	asly_args.add_argument("--spades-k-list",
		help="SPAdes: List of kmers, comma-separated even numbers e.g. '21,33,55,77'")
	# Unicycler options
	asly_args.add_argument("--unicycler-mode", default="normal",
		choices=['conservative', 'normal', 'bold'],
		help="Unicycler: assember mode: conservative, normal (default) or bold.")
	asly_args.add_argument("--linear-seqs", default=0, type=int,
		help="Expected number of linear sequences")
	asly_args.add_argument("--extract-replicons", action="store_true",
		help="Unicycler: extract complete replicons (e.g. plasmids)"
		+ " from the short-read based assembly to separate files")
	# Flye options
	asly_args.add_argument("--flye-short-polish", action="store_true",
		help="Perform polishing of Flye assembly with short reads using racon.")
	asly_args.add_argument("--flye-skip-long-polish", action="store_true",
		help="Skip stage of genome polishing with long reads.")
	asly_args.add_argument("--perform-polishing", action="store_true",
		help="Perform polishing. Useful only for flye assembly of long reads and short reads available.")
	asly_args.add_argument("--polishing-iterations", default=1, type=int,
		help="Number of polishing iterations, default: 1.")

	check_args = parser.add_argument_group(title="Genome check settings")
	# phiX
	check_args.add_argument("--check-phix", action="store_true",
		help="Check genome for presence of PhiX control sequence.")
	# CheckM
	check_args.add_argument("--checkm-mode",
		default="taxonomy_wf", choices=['taxonomy_wf', 'lineage_wf'],
		help="CheckM working mode. Default is checking for domain-specific marker-set.")
	check_args.add_argument("--checkm-rank",
		help="Rank of taxon for CheckM. Run 'checkm taxon_list' for details.")
	check_args.add_argument("--checkm-taxon",
		help="Taxon for CheckM. Run 'checkm taxon_list' for details.")
	check_args.add_argument("--checkm-full-tree", action="store_true",
		help="Use full tree for inference of marker set, requires LOTS of memory.")

	anno_args = parser.add_argument_group(title="Annotation settings")
	anno_args.add_argument("-g", "--genome",
		help="Genome assembly (when starting with assembled genome).")
	anno_args.add_argument("--locus-tag",
		help="Locus tag prefix. If not provided prefix will be by bakta.")
	anno_args.add_argument("--compliant", action='store_true',
		help="Force Genbank/ENA/DDJB compliance.")
	anno_args.add_argument("--minimum-contig-length", type=int,
		help="Minimum sequence length in genome assembly.")
	anno_args.add_argument("--prefix",
		help="Prefix for annotated files.")
	anno_args.add_argument("--translation-table", default=11, type=int,
		choices=[4, 11], help="Translation table: 11/4, default: 11.")

	args = parser.parse_args()

	if args.pe_1 and args.pe_2:
		assert args.pe_1 != args.pe_2, 'Same file(s) for forward and reverse reads!'

	if args.mp_1 and args.mp_2:
		assert args.mp_1 != args.mp_2, 'Same file(s) for forward and reverse mate-pair reads!'

	if (args.assembler == 'spades'
		and not (
			bool(args.pe_1)
			or bool(args.pe_merged)
			or bool(args.single_end)
			or bool(args.mp_1)
		)
	):
		logger.error("Impossible to run SPAdes without short reads!")
		raise Exception("Insufficient input data.")

	if args.assembler == 'flye':

		if args.nanopore and args.pacbio:
			logger.error("Impossible to run Flye on mixed long reads!")
			raise Exception("Incorrect input data.")

		if not (bool(args.nanopore) or bool(args.pacbio)):
			logger.error("Impossible to run Flye without long reads!")
			raise Exception("Insufficient input data.")

		if not args.genome_size_estimation:
			args.calculate_genome_size = True
			logger.info("Genome size was not provided. It will be calculated with mash.")

	if (args.flye_short_polish
		or args.first_step == 'polishing'
		or args.last_step == 'polishing'
	):
		args.perform_polishing = True

	args.output_dir = os.path.abspath(args.output_dir)

	return args


def check_reads(args):
	'''Check reads provided by user.

	Returns:
	reads (list(dict)) : paths to all existing reads.
	'''
	libraries = []
	read_list = [args.pe_1, args.pe_2, args.single_end, args.pe_merged,
		args.mp_1, args.mp_2, args.pacbio, args.nanopore]
	read_list = list(filter(None.__ne__, read_list))

	logger.info("Checking input files.")

	for reads in read_list:
		for f_name in reads:
			if f_name != "n/a" and not os.path.isfile(f_name):
				logger.error("File %s doesn't exist", f_name)
				raise FileNotFoundError(f"File {f_name} doesn't exist")

	short_libs = {}
	if args.pe_1 and args.pe_2:
		short_libs["forward"] = args.pe_1
		short_libs["reverse"] = args.pe_2
	if args.single_end:
		short_libs["single"] = args.single_end
	if args.pe_merged:
		short_libs["merged"] = args.pe_merged

	if len(short_libs.values()) > 0:
		short_lib_count = max(list(map(len, short_libs.values())))
		for i in range(short_lib_count):
			libraries.append({"type": "short"})
			for lib_type, lib in short_libs.items():
				if len(lib) > i and lib[i] != "n/a":
					libraries[i][lib_type] = os.path.abspath(lib[i])

	if args.mp_1 and args.mp_2:
		for pair in list(zip(args.mp_1, args.mp_2)):
			libraries.append({
				"type": "mate-pair",
				"forward": os.path.abspath(pair[0]),
				"reverse": os.path.abspath(pair[1])
			})

	if args.pacbio:
		for lib in args.pacbio:
			libraries.append({
				"type": "pacbio",
				"single": os.path.abspath(lib)
			})

	if args.nanopore:
		for lib in args.nanopore:
			libraries.append({
				"type": "nanopore",
				"single": os.path.abspath(lib)
			})

	if len(libraries) == 0:
		logger.critical("No reads provided for genome assembly")
		raise Exception("No reads provided for genome assembly")

	return libraries


def read_qc(args, reads):
	'''Perform read QC with fastp'''
	logger.info("Read quality control started")
	qcoutdir = create_subdir(args.output_dir, "readQC")
	precmd = ["fastp", "-L", "-Q", "-G", "-A", "-z", "1", "--stdout",
		"-w", str(args.threads)]
	for lib in reads:
		for key, value in lib.items():
			if key == "type":
				continue
			prefix = os.path.join(qcoutdir, os.path.split(value)[-1])
			cmd = precmd + ["-i", value, "-h", f"{prefix}.html", "-j", f"{prefix}.json"]
			logger.debug("QC of %s", value)
			run_external(args, cmd)


def mash_estimate(args, reads):
	'''Estimates genome size using mash.

	Parameters:
	args (argparse.Namespace)
	reads (dict) : input reads

	Returns:
	estimated genome size (int) or None

	Parsing mash output to extract estimated genome size from stderr lines:

	Estimated genome size: 1.234e+06
	Estimated coverage:    56.789

	Genome size with highest coverage is taken from the sorted list of tuples (size, coverage).
	'''
	reads_to_sketch = []

	for library in reads:
		if library['type'] == "mate-pair":
			continue
		for readfile in [v for k, v in library.items() if k != "type"]:
			reads_to_sketch.append(readfile)

	if len(reads_to_sketch) == 0:
		logger.error("Not possible to estimate gemome size: reads missing")
		return None

	sketchprefix = os.path.join(os.path.dirname(reads_to_sketch[0]), "sketch")
	cmd = ["mash", "sketch", "-r", "-m", str(args.mash_kmer_copies), "-o", sketchprefix]
	cmd += reads_to_sketch
	logger.info("Estimating genome size with mash using: %s",
		", ".join(reads_to_sketch))
	r = run_external(args, cmd, keep_stderr=True)

	if r:
		result = r.stderr.split("\n")
		estimations = [float(x.split()[-1]) for x in result if "Estimated" in x.split()]
		best_estimation = sorted(zip(estimations[::2], estimations[1::2]), key=lambda x: -x[1])[0]
		logger.info("Estimated genome size is %s bp at coverage %s.",
			int(best_estimation[0]), best_estimation[1])
		return int(best_estimation[0])
	logger.error("Genome size estimation with \"mash\" failed.")
	return None


def map_short_reads(args, assembly, reads, target):
	'''Short read mapping with minimap2.

	Parameters:
	args
	reads (str) : path to reads (single file) in FASTQ format [gzippd]
	assembly (str) : path to assembly in FASTA format
	target (str) : path to file to store mapping

	Returns:
	target (str) : Path to mapping file
	'''
	cmd = ["minimap2",
		"-x", "sr", "-t", str(args.threads),
		"-a", "-o", target, assembly, reads]

	logger.info("Mapping reads: %s", reads)
	if run_external(args, cmd):
		return target
	logger.error("Unsuccesful mapping.")
	return None


def racon_polish(args, assembly, reads) -> str:
	'''Genome assembly polishing with racon

	Returns:
	assembly (str) : path to polished assembly
	'''
	polish_dir = os.path.join(args.output_dir, "polishing")
	if not os.path.isdir(polish_dir):
		polish_dir = create_subdir(args.output_dir, "polishing")
	target = os.path.join(polish_dir, "mapping.sam")
	for library in [x for x in reads if x['type'] == "short"]:
		for readfile in [x for x in library.values() if x != "short"]:
			logger.debug("Racon genome polishing with: %s", readfile)
			mapping = map_short_reads(args, assembly, readfile, target)
			if mapping:
				cmd = ["racon", "-t", str(args.threads), readfile, mapping, assembly]
				r = run_external(args, cmd, keep_stdout=True)
				if os.path.exists(mapping):
					os.remove(mapping)
				if r:
					digest = hashlib.md5(r.stdout.encode('utf-8')).hexdigest()
					suffix = "".join(
						[chr(65 + (int(digest[x], 16) + int(digest[x + 1], 16)) % 26) for x in range(0, 20, 2)]
					)
					fname = os.path.join(polish_dir, f"polished.{suffix}.fna")
					try:
						with open(fname, 'w') as handle:
							handle.write(r.stdout)
							assembly = fname
					except Exception:
						logger.error("Error during polishing: impossible to write file %s", fname)
			else:
				logger.error("Impossible to perform polishing.")

	return assembly


def annotate(args) -> str:
	'''Returns path to annotated genome (FASTA)'''
	logger.info("Genome annotation started")

	try:
		version_stdout = subprocess.run(
			["bakta", "--version"], encoding="utf-8", check=True,
			stderr=subprocess.PIPE, stdout=subprocess.PIPE).stdout
		version = re.search(r'bakta (\S)', version_stdout)[1]
		logger.debug("bakta version %s available.", version)
	except Exception as e:
		logger.critical("Failed to run bakta!")
		raise e

	annodir = os.path.join(args.output_dir, "annotation")

	cmd = ["bakta", args.genome, "-o", annodir,
		'--genus', args.genus, '--species', args.species,
		"--threads", str(args.threads),
		'--translation-table', str(args.translation_table)]
	if args.prefix:
		cmd += ['--prefix', args.prefix]
	if args.strain:
		cmd += ['--strain', args.strain]
	if args.locus_tag:
		cmd += ['--locus-tag', args.locus_tag]
	if args.compliant:
		cmd += ['--compliant']
	if args.minimum_contig_length:
		cmd += ['--min-contig-length', str(args.minimum_contig_length)]

	if run_external(args, cmd):
		stem = Path(args.genome).stem
		args.genome = os.path.join(annodir, f"{stem}.fna")
	return args.genome


def check_phix(args):
	'''Checks assembly for phiX sequence

	Returns
	Path to clean genome assembly
	'''

	logger.info("Checking assembly for presence of Illumina phiX control")
	phix_path = os.path.join(
		os.path.dirname(os.path.abspath(__file__)),
		"data/phiX174.fasta"
	)
	blast_format = "6 sseqid pident slen length"
	cmd = ["blastn",
		"-query", phix_path, "-subject", args.genome,
		"-outfmt", blast_format, "-evalue", "1e-6"]

	logger.debug("Running: %s", " ".join(cmd))
	try:
		blast_out = str(subprocess.run(
			cmd, universal_newlines=True, check=True,
			stderr=subprocess.PIPE, stdout=subprocess.PIPE
		).stdout).rstrip()
	except Exception as e:
		logger.error("Error during BLAST+ run")
		raise e

	phix_contigs = []
	if bool(blast_out):
		for line in blast_out.split('\n'):
			sseqid, pident, slen, length = line.split("\t")
			if float(pident) > 95.0 and int(slen) / int(length) > 0.5:
				phix_contigs.append(sseqid)
		phix_contigs = list(set(phix_contigs))

		if len(phix_contigs) > 0:
			logger.info("PhiX was found in: %s", ', '.join(phix_contigs))
			newgenome = os.path.join(args.output_dir, "assembly.nophix.fasta")
			records = [x for x in SeqIO.parse(args.genome, "fasta") if x.id not in phix_contigs]
			with open(newgenome, "w") as handle:
				SeqIO.write(records, handle, "fasta")
				args.genome = newgenome
		else:
			logger.info("PhiX wasn't found.")

	return args.genome


def run_checkm(args):
	'''Run CheckM

	Returns
	Path to CheckM output file (str) or None
	'''
	checkm_indir = create_subdir(args.output_dir, "checkm_tmp_in")
	try:
		shutil.copy(args.genome, checkm_indir)
	except Exception as e:
		logger.critical("Failed to copy \"%s\" to %s", args.genome, checkm_indir)
		raise e

	checkm_outdir = os.path.join(args.output_dir, "checkm")
	checkm_outfile = os.path.join(args.output_dir, "checkm.result.txt")
	checkm_ext = os.path.splitext(args.genome)[1]

	if args.checkm_mode == "taxonomy_wf":

		cmd = ["checkm", "taxonomy_wf",
			"-f", checkm_outfile,
			"-x", checkm_ext,
			"-t", str(args.threads)]

		try:
			checkm_taxon_list = subprocess.run(
				["checkm", "taxon_list"],
				encoding="utf-8",
				universal_newlines=True,
				stderr=subprocess.DEVNULL,
				stdout=subprocess.PIPE,
				check=True
			).stdout
		except Exception as e:
			logger.critical("Failed to run CheckM!")
			raise e

		if args.checkm_taxon and args.checkm_rank:
			found = False
			for line in checkm_taxon_list.split("\n"):
				if args.checkm_taxon in line and args.checkm_rank in line:
					found = True
					break
			if not found:
				logger.error("Taxon %s of rank %s not available for CheckM",
					args.checkm_taxon, args.checkm_rank)
				args.checkm_taxon = None
				args.checkm_rank = None

		if not args.checkm_taxon or not args.checkm_rank:
			args.checkm_taxon = args.domain
			args.checkm_rank = "domain"

		logger.info("%s marker set will be used for CheckM", args.checkm_taxon)

		cmd += [args.checkm_rank, args.checkm_taxon, checkm_indir, checkm_outdir]

	else:
		cmd = ["checkm", "lineage_wf",
			"-f", checkm_outfile,
			"-t", str(args.threads),
			"-x", checkm_ext]
		if not args.checkm_full_tree:
			cmd += ["--reduced_tree"]
		cmd += ["--pplacer_threads", str(args.threads),
			checkm_indir, checkm_outdir]

	r = run_external(args, cmd)

	# Cleaning after CheckM
	shutil.rmtree(checkm_indir)
	shutil.rmtree(checkm_outdir)

	if r:
		return checkm_outfile
	logger.error("CheckM didn't finish properly.")
	return None


def check_last_step(args, step):
	if args.last_step == step:
		logger.info("Workflow finished!")
		sys.exit(0)


def create_out_dir(args):
	if os.path.isdir(args.output_dir):
		if args.force:
			shutil.rmtree(args.output_dir)
		else:
			logger.critical(
				"Output directory \"%s\" already exists. Use --force to overwrite.",
				args.output_dir)
			raise FileExistsError(args.output_dir)
	try:
		os.mkdir(args.output_dir)
	except Exception as e:
		logger.critical("Imposible to create directory \"%s\"!", args.output_dir)
		raise e
	return os.path.abspath(args.output_dir)


def main():
	global logger
	logger = logging.getLogger("main")
	logger.setLevel(logging.DEBUG)
	formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
	ch = logging.StreamHandler()
	ch.setLevel(logging.INFO)
	ch.setFormatter(formatter)
	logger.addHandler(ch)

	args = parse_args()

	args.output_dir = create_out_dir(args)
	logger.debug("Writing output to %s", create_out_dir(args))

	zgalogfile = os.path.join(args.output_dir, "zga.log")
	fh = logging.FileHandler(zgalogfile)
	fh.setLevel(logging.DEBUG)
	fh.setFormatter(formatter)
	logger.addHandler(fh)

	logger.info("ZGA ver. %s", __version__)
	logger.info("Full log location: %s", zgalogfile)

	steps = {
		"readqc": 1,
		"processing": 2,
		"assembling": 3,
		"polishing": 4,
		"check_genome": 5,
		"annotation": 6
	}
	args.first_step = steps[args.first_step]
	args.last_step = steps[args.last_step]
	if args.first_step <= 4:
		reads = check_reads(args)
		logger.debug("Input reads: %s", reads)

	if (args.first_step > 3
		and (
			not args.genome
			or not os.path.isfile(args.genome)
		)
	):
		logger.error("Genome assembly is not provided")
		raise FileNotFoundError()

	# Read QC
	if args.first_step == 1:
		read_qc(args, reads)
		check_last_step(args, 1)

	# Processing
	if args.first_step <= 2:
		reads = read_processing(args, reads)
		logger.debug("Processed reads: %s", str(reads))
		check_last_step(args, 2)

	# Assembly
	if args.first_step <= 3:
		if args.calculate_genome_size:
			estimated_genome_size = mash_estimate(args, reads)
		else:
			estimated_genome_size = args.genome_size_estimation
		args.genome = assemble(args, reads, estimated_genome_size)
		check_last_step(args, 3)

	# Short read polishing is only meaningful for flye assembly
	if args.first_step <= 4 and args.perform_polishing:
		short_reads_list = [lib for lib in reads if lib['type'] == "short"]
		if len(short_reads_list) > 0:
			for x in range(args.polishing_iterations):
				logger.info("Performing genome polishing. Iteration %s.", x)
				args.genome = racon_polish(args, args.genome, short_reads_list)
			args.genome = shutil.copy2(args.genome,
				os.path.join(args.output_dir, "assembly.polished.fasta"))
			logger.info("Genome polishing finished. Polished genome: %s.", args.genome)
		else:
			logger.error("Short reads are not provided for polishing of flye assembly.")
		check_last_step(args, 4)

	# Assembly QC
	if args.first_step <= 5:
		logger.info("Checking genome quality")
		if args.check_phix:
			args.genome = check_phix(args)
		checkm_outfile = run_checkm(args)
		if checkm_outfile:
			with open(checkm_outfile) as result:
				completeness, contamination, heterogeneity = list(
					map(float, result.readlines()[3].split()[-3::1])
				)
				if completeness < 80.0 or (completeness - 5.0 * contamination < 50.0):
					logger.warning("The genome assembly has low quality!")
				logger.info("Genome completeness: %s%%", completeness)
				logger.info("Genome contamination: %s%%", contamination)
				logger.info("Genome heterogeneity: %s%%", heterogeneity)
		check_last_step(args, 5)

	# Annotation
	if args.first_step <= 6:
		if not args.genome or not os.path.isfile(args.genome):
			logger.error("Genome assembly is not provided")
			raise FileNotFoundError()
		args.genome = annotate(args)
		check_last_step(args, 6)


if __name__ == "__main__":
	main()
