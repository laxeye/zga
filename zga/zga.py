#!/usr/bin/env python3
import argparse
import os.path
import logging
import sys
import shutil
import subprocess
import re
from Bio import SeqIO
import hashlib


def parse_args():
	parser = argparse.ArgumentParser(description="ZGA genome assembly and annotation pipeline")

	# General options
	general_args = parser.add_argument_group(title="General options", description="")
	general_args.add_argument("-s", "--first-step", help="First step of the pipeline", default="readqc",
		choices=["readqc", "processing", "assembling", "polishing", "check_genome", "annotation"])
	general_args.add_argument("-l", "--last-step", help="Last step of the pipeline", default="annotation",
		choices=["readqc", "processing", "assembling", "polishing", "check_genome", "annotation"])
	general_args.add_argument("-o", "--output-dir", required=True, help="Output directory")
	general_args.add_argument("--force", action="store_true",
		help="Overwrite output directory if exists")
	# parser.add_argument("--tmp-dir", default="zga-temp", help="Temporary directory")
	general_args.add_argument("-t", "--threads", type=int, default=1,
		help="Number of CPU threads to use (where possible)")
	general_args.add_argument("-m", "--memory-limit", type=int, default=8,
		help="Memory limit in GB (default 8)")
	general_args.add_argument("--genus", default="Unknown", help="Provide genus if known")
	general_args.add_argument("--species", default="sp.", help="Provide species if known")
	general_args.add_argument("--strain", help="Provide strain if known")
	general_args.add_argument("--transparent", action="store_true",
		help="Show output from tools inside the pipeline")
	general_args.add_argument("--domain", default="bacteria", choices=['archaea', 'bacteria'],
		help="Provide prokaryotic domain: bacteria or archaea")

	# Input
	input_args = parser.add_argument_group(title="Input files and options",
		description="Sequencing reads should be in FASTQ format and may be GZipped.")
	input_args.add_argument("-1", "--pe-1", help="FASTQ file with first (left) paired-end reads")
	input_args.add_argument("-2", "--pe-2", help="FASTQ file with second (right) paired-end reads")
	input_args.add_argument("--pe-merged", help="FASTQ file  with merged overlapped paired-end reads")
	input_args.add_argument("-S", "--single-end", help="FASTQ file with unpaired or single-end reads")
	input_args.add_argument("--mp-1", help="Mate pair forward reads. SPAdes only")
	input_args.add_argument("--mp-2", help="Mate pair forward reads. SPAdes only")
	# parser.add_argument("--pe-interleaved", help="Pair-end interleaved reads")
	# parser.add_argument("--mp-interleaved", help="Mate pair interleaved reads")
	input_args.add_argument("--pacbio", help="PacBio reads")
	input_args.add_argument("--nanopore", help="Nanopore reads")

	# Read processing
	reads_args = parser.add_argument_group(title="Read processing settings")
	reads_args.add_argument("-q", "--quality-cutoff", type=int, default=25,
		help="Base quality cutoff for short reads, default: 25")
	reads_args.add_argument("--adapters",
		help="Adapter sequences for short reads trimming (FASTA). "
		+ "By default Illumina adapter sequences are used.")
	reads_args.add_argument("--filter-by-tile", action="store_true",
		help="Filter Illumina reads based on positional quality over a flowcell.")
	reads_args.add_argument("--min-short-read-length", type=int, default=33,
		help="Minimum short read length to keep after quality trimming.")
	reads_args.add_argument("--bbmerge-extend", type=int,
		help="Perform read extension by specified length based on k-mers "
		+ "if initial merging wasn't succesfull.")
	reads_args.add_argument("--bbmerge-extend-kmer", type=int, default=40,
		help="K-mer length for read extension, default 40.")
	reads_args.add_argument("--bbmerge-trim", type=int,
		help="Before merging trim bases with phred score less than a specified value.")
	reads_args.add_argument("--genome-size-estimation", action="store_true",
		help="Estimate genome size with mash.")
	reads_args.add_argument("--estimated-genome-size",
		help="Estimated genome for Flye assembler, if known.")
	reads_args.add_argument("--mash-kmer-copies", type=int, default=10,
		help="Minimum copies of each k-mer to inslude in size estimation")
	# Mate pair read processing
	reads_args.add_argument("--use-unknown-mp", action="store_true",
		help="Include reads that are probably mate pairs (default: only known MP used)")

	# Assembly
	asly_args = parser.add_argument_group(title="Assembly settings")
	asly_args.add_argument("-a", "--assembler", default="unicycler", choices=["spades", "unicycler", "flye"],
		help="Assembler: unicycler (default; better quality), spades (faster, may use mate-pair reads) or Flye (long reads only).")
	asly_args.add_argument("--no-correction", action="store_true",
		help="Disable read correction in SPAdes")
	# Spades options
	asly_args.add_argument("--use-scaffolds", action="store_true",
		help="SPAdes: Use assembled scaffolds.")
	asly_args.add_argument("--spades-k-list",
		help="SPAdes: List of kmers, comma-separated even numbers e.g. '21,33,55,77'")
	# Unicycler options
	asly_args.add_argument("--unicycler-mode", default="normal",
		choices=['conservative', 'normal', 'bold'],
		help="Unicycler: assember mode: conservative, normal (default) or bold.")
	asly_args.add_argument("--linear-seqs", default=0,
		help="Expected number of linear sequences")
	asly_args.add_argument("--extract-replicons", action="store_true",
		help="Unicycler: extract replicons (e.g. plasmids) from the short-read based assembly to separate files")
	# Flye options
	asly_args.add_argument("--flye-short-polish", action="store_true",
		help="Perform polishing of Flye assembly with short reads using racon.")
	asly_args.add_argument("--skip-flye-long-polish", action="store_true",
		help="Skip stage of genome polishing with long reads.")
	asly_args.add_argument("--perform-polishing", action="store_true",
		help="Perform polishing. Useful only for flye assembly of long reads and short reads available.")
	asly_args.add_argument("--polishing-iterations", default=1,
		help="Number of polishing iterations.")

	check_args = parser.add_argument_group(title="Genome check settings")
	# phiX
	check_args.add_argument("--check-phix", action="store_true",
		help="Check genome for presence of PhiX control sequence.")
	# CheckM
	check_args.add_argument("--checkm-mode", default="taxonomy_wf", choices=['taxonomy_wf', 'lineage_wf'],
		help="Select CheckM working mode. Default is checking for domain-specific marker-set.")
	check_args.add_argument("--checkm-rank", help="Rank of taxon for CheckM. Run 'checkm taxon_list' for details.")
	check_args.add_argument("--checkm-taxon", help="Taxon for CheckM. Run 'checkm taxon_list' for details.")
	check_args.add_argument("--checkm-full-tree", action="store_true",
		help="Use full tree for inference of marker set, requires LOTS of memory.")

	anno_args = parser.add_argument_group(title="Annotation settings")
	# Annotation
	anno_args.add_argument("-g", "--genome", help="Genome assembly (when starting from annotation).")
	anno_args.add_argument("--gcode", default=11, type=int, help="Genetic code.")
	anno_args.add_argument("--locus-tag",
		help="Locus tag prefix. If not provided prefix will be generated from MD5 checksum.")
	anno_args.add_argument("--locus-tag-inc", default=10, type=int,
		help="Locus tag increment, default = 10")
	anno_args.add_argument("--center-name", help="Genome center name.")
	anno_args.add_argument("--minimum-contig-length", help="Minimum sequence length in genome assembly.")

	args = parser.parse_args()

	if args.assembler == 'spades' and not ( bool(args.pe_1) |
		bool(args.pe_merged) | bool(args.single_end) | bool(args.mp_1) ) :
		logger.error("Impossible to run SPAdes without short reads!")
		raise Exception("Bad parameters.")

	if args.assembler == 'flye':

		if args.nanopore and args.pacbio:
			logger.error("Impossible to run Flye on mixed long reads!")
			raise Exception("Bad parameters.")

		if not ( bool(args.nanopore) | bool(args.pacbio) ):
			logger.error("Impossible to run Flye without long reads!")
			raise Exception("Bad parameters.")

		if not args.estimated_genome_size:
			args.genome_size_estimation = True
			logger.info("Genome size was not provided. It will be calculated with mash.")

	if args.flye_short_polish:
		args.perform_polishing = True

	return args


def check_reads(args):
	reads = {}
	read_list = [args.pe_1, args.pe_2, args.single_end, args.pe_merged,
	args.mp_1, args.mp_2, args.pacbio, args.nanopore]

	logger.info("Checking input files.")

	for f in read_list:
		if f and not os.path.isfile(f):
			logger.error("File %s doesn't exist", f)
			raise FileNotFoundError("File %s doesn't exist" % f)

	if args.pe_1 and args.pe_2:
		reads['pe'] = (os.path.abspath(args.pe_1), os.path.abspath(args.pe_2))
	elif args.pe_1 or args.pe_2:
		logger.error("Single end reads provided as paired. Please use \"--single-end\" option")
		sys.exit(1)

	if args.mp_1 and args.mp_2:
		reads['mp'] = (os.path.abspath(args.mp_1), os.path.abspath(args.mp_2))

	if args.pe_merged:
		reads['merged'] = os.path.abspath(args.pe_merged)

	if args.single_end:
		reads['single'] = os.path.abspath(args.single_end)

	if args.pacbio:
		reads['pacbio'] = os.path.abspath(args.pacbio)

	if args.nanopore:
		reads['nanopore'] = os.path.abspath(args.nanopore)

	return reads


def create_subdir(parent, child):
	path = os.path.join(parent, child)
	try:
		os.mkdir(path)
	except Exception as e:
		logger.critical("Impossible to create directory \"%s\"", path)
		raise e
	return path


def run_external(args, cmd, keep_stdout=False, keep_stderr=False):
	logger.debug("Running: %s", " ".join(cmd))
	stderr_dest = subprocess.PIPE if keep_stderr else subprocess.DEVNULL
	stdout_dest = subprocess.PIPE if keep_stdout else subprocess.DEVNULL

	if args.transparent:
		if keep_stdout:
			r = subprocess.run(cmd, stdout=subprocess.PIPE, encoding="utf-8")
		elif keep_stderr:
			r = subprocess.run(cmd, stderr=subprocess.PIPE, encoding="utf-8")
		else:
			r = subprocess.run(cmd, stdout=subprocess.DEVNULL)
	else:
		r = subprocess.run(cmd, stderr=stderr_dest, stdout=stdout_dest, encoding="utf-8")

	if r.returncode != 0:
		logger.error("Non-zero return code of %s", " ".join(cmd))

	return r


def read_qc(args, reads):
	logger.info("Read quality control started")
	qcoutdir = create_subdir(args.output_dir, "readQC")
	precmd = ["fastp", "-L", "-Q", "-G", "-A", "-z", "1", "--stdout", "-w", str(args.threads)]
	for r in reads.values():
		if isinstance(r, (list, tuple)):
			prefix = os.path.join(qcoutdir,os.path.split(r[0])[-1])
			cmd = precmd + ["-i", r[0], "-I", r[1]]
		else:
			prefix = os.path.join(qcoutdir,os.path.split(r)[-1])
			cmd = precmd + ["-i", r]
		cmd += ["-h", f"{prefix}.html", "-j", f"{prefix}.json"]
		logger.debug("QC of %s", r)
		run_external(args, cmd)


def filter_by_tile(args, reads, readdir):
	filtered_pe_r1 = os.path.join(readdir, "filtered_pe_r1.fq.gz")
	filtered_pe_r2 = os.path.join(readdir, "filtered_pe_r2.fq.gz")

	cmd = ["filterbytile.sh", f"in={reads['pe'][0]}", f"in2={reads['pe'][1]}",
	f"out={filtered_pe_r1}", f"out2={filtered_pe_r2}"]

	if run_external(args, cmd).returncode == 0:
		for f in reads['pe']:
			if os.path.dirname(f) == readdir and os.path.exists(f):
				os.remove(f)
		reads['pe'] = (filtered_pe_r1, filtered_pe_r2)

	return reads


def merge_bb(args, reads, readdir):
	'''Performs overlapping paired reads merging (and extension).
	reads - dict of input reads
	readdir - path to reads output directory
	Retruns reads, modified if bbmerge returns 0.
	'''
	notmerged_r1 = os.path.join(readdir, "nm.pe_1.fq.gz")
	notmerged_r2 = os.path.join(readdir, "nm.pe_2.fq.gz")
	merged = os.path.join(readdir, "merged.fq.gz")
	cmd = ["bbmerge.sh", f"Xmx={args.memory_limit}G", f"t={args.threads}",
		f"in1={reads['pe'][0]}", f"in2={reads['pe'][1]}",
		f"outu1={notmerged_r1}", f"outu2={notmerged_r2}", f"out={merged}"]
	if args.bbmerge_extend:
		cmd += [f"extend2={args.bbmerge_extend}",
		f"k={args.bbmerge_extend_kmer}", "rsem=t"]
	if args.bbmerge_trim:
		cmd += ["qtrim2=t", f"trimq={args.bbmerge_trim}"]
	logger.info("Merging paired-end reads.")

	if run_external(args, cmd).returncode == 0:
		for f in reads['pe']:
			if os.path.dirname(f) == readdir and os.path.exists(f):
				os.remove(f)
		reads['merged'] = merged
		reads['pe'] = (notmerged_r1, notmerged_r2)

	return reads


def bbduk_process(args, reads, readdir):
	bbduk_kmer = 19  # K-mer length for contaminant/adapter removal
	precmd = ["bbduk.sh", f"Xmx={args.memory_limit}G", f"t={args.threads}",
			f"ref={args.adapters}", f"k={bbduk_kmer}", "ktrim=r",
			"qtrim=r", f"trimq={args.quality_cutoff}", "entropy=0.1",
			f"minlength={args.min_short_read_length}"]

	if "pe" in reads.keys():
		logger.info("Trimming and filtering paired end reads")
		out_pe1 = os.path.join(readdir, "pe_1.fq")
		out_pe2 = os.path.join(readdir, "pe_2.fq")
		out_stats = os.path.join(readdir, "bbduk.stats.pe.txt")
		cmd = precmd + [f"in={reads['pe'][0]}", f"in2={reads['pe'][1]}",
			f"out={out_pe1}", f"out2={out_pe2}", f"stats={out_stats}"]

		if run_external(args, cmd).returncode == 0:
			for f in reads['pe']:
				if os.path.dirname(f) == readdir and os.path.exists(f):
					os.remove(f)
			reads['pe'] = (out_pe1, out_pe2)

	if "single" in reads.keys():
		logger.info("Trimming and filtering single end reads")
		out = os.path.join(readdir, "single.fq")
		out_stats = os.path.join(readdir, "bbduk.stats.se.txt")
		cmd = precmd + [f"in={reads['single']}", f"out={out}",
			f"stats={out_stats}"]

		if run_external(args, cmd).returncode == 0:
			reads['single'] = out

	if "merged" in reads.keys():
		logger.info("Trimming and filtering merged paired-end reads")
		out = os.path.join(readdir, "merged.fq")
		out_stats = os.path.join(readdir, "bbduk.stats.merged.txt")
		cmd = precmd + [f"in={reads['merged']}", f"out={out}",
			f"stats={out_stats}"]

		if run_external(args, cmd).returncode == 0:
			reads['merged'] = out

	return reads


def read_processing(args, reads):
	logger.info("Reads processing started")
	readdir = create_subdir(args.output_dir, "reads")
	illumina_adapters = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data/illumina.adapters.fasta")

	if args.adapters and os.path.isfile(args.adapters):
		args.adapters = os.path.abspath(args.adapters)
	else:
		args.adapters = illumina_adapters

	if args.filter_by_tile and "pe" in reads.keys():
		reads = filter_by_tile(args, reads, readdir)

	reads = bbduk_process(args, reads, readdir)

	# Merging overlapping paired-end reads
	if "merged" not in reads.keys() and "pe" in reads.keys():
		reads = merge_bb(args, reads, readdir)

	# Processing Illumina mate-pairs
	if "mp" in reads.keys():
		reads = mp_read_processing(args, reads, readdir)

	logger.info("Read processing finished")

	return reads


def mash_estimate(args, reads):
	reads_to_sketch = []

	for readfile in reads.values():
		if isinstance( readfile, (list, tuple) ):
			reads_to_sketch += list(readfile)
		else:
			reads_to_sketch.append(readfile)

	if len(reads_to_sketch) == 0:
		logger.error("Not possible to estimate gemome size: reads missing")
		return None

	sketchprefix = os.path.join(os.path.dirname(reads_to_sketch[0]), "sketch")
	cmd = ["mash", "sketch", "-r", "-m", str(args.mash_kmer_copies), "-o", sketchprefix]
	cmd += reads_to_sketch
	logger.info("Estimating genome size with mash using: %s", ", ".join(reads_to_sketch))
	r = run_external(args, cmd, keep_stderr=True)

	'''
	Parsing mash output to extract estimated genome size
	from lines looking like:
	Estimated genome size: 1.234e+06
	Estimated coverage:    23.644
	We create list of tuples (size, coverage), sort it and get value with greatest covarage.
	'''
	if r.returncode == 0:
		result = r.stderr.split("\n")
		estimations = [float(x.split()[-1]) for x in result if "Estimated" in x.split()]
		best_estimation = sorted(zip(estimations[::2], estimations[1::2]), key=lambda x: -x[1])[0]
		logger.info("Estimated genome size is %s bp at coverage %s.",
			int(best_estimation[0]), best_estimation[1])
		return int(best_estimation[0])
	else:
		logger.error("Genome size estimation with \"mash\" failed.")
		return None


def mp_read_processing(args, reads, readdir):
	prefix = os.path.join(readdir, "nxtrim")

	cmd = ["nxtrim", "-1", reads['mp'][0], "-2", reads['mp'][1], "--separate"]
	cmd += ["--justmp", "-O", prefix, "-l", str(args.min_short_read_length)]
	logger.info("Processing mate-pair reads.")
	if run_external(args, cmd).returncode == 0:
		if args.use_unknown_mp:
			with open(f"{prefix}_R1.all.fastq.gz", "wb") as dest:
				with open(f"{prefix}_R1.mp.fastq.gz", "rb") as src:
					shutil.copyfileobj(src, dest)
				with open(f"{prefix}_R1.unknown.fastq.gz", "rb") as src:
					shutil.copyfileobj(src, dest)
			with open(f"{prefix}_R2.all.fastq.gz", "wb") as dest:
				with open(f"{prefix}_R2.mp.fastq.gz", "rb") as src:
					shutil.copyfileobj(src, dest)
				with open(f"{prefix}_R2.unknown.fastq.gz", "rb") as src:
					shutil.copyfileobj(src, dest)
			reads['mp'] = (f"{prefix}_R1.all.fastq.gz", f"{prefix}_R2.all.fastq.gz")
		else:
			reads['mp'] = (f"{prefix}_R1.mp.fastq.gz", f"{prefix}_R2.mp.fastq.gz")

	return reads


def map_short_reads(args, assembly, reads, target):
	cmd = ["minimap2", "-x", "sr", "-t", str(args.threads), "-a", "-o", target, assembly, reads]

	logger.info("Mapping reads: %s", reads)
	if run_external(args, cmd).returncode == 0:
		return target
	else:
		logger.error("Unsuccesful mapping.")
		return None


def racon_polish(args, assembly, reads):
	polish_dir = os.path.join(args.output_dir, "polishing")
	if not os.path.isdir(polish_dir):
		polish_dir = create_subdir(args.output_dir, "polishing")
	target = os.path.join(polish_dir, "mapping.sam")
	for readfile in reads:
		if isinstance(readfile, (list, tuple)):
			assembly = racon_polish(args, assembly, readfile)
		else:
			logger.debug("Racon genome polishing with: %s", readfile)
			mapping = map_short_reads(args, assembly, readfile, target)
			if mapping:
				cmd = ["racon", "-t", str(args.threads), readfile, mapping, assembly]
				r = run_external(args, cmd, keep_stdout=True)
				if os.path.exists(mapping):
					os.remove(mapping)
				if r.returncode == 0:
					digest = hashlib.md5(r.stdout.encode('utf-8')).hexdigest()
					suffix = "".join([chr(65 + (int(digest[x], 16) + int(digest[x + 1], 16)) % 26) for x in range(0, 20, 2)])
					fname = os.path.join(polish_dir, f"polished.{suffix}.fna")
					try:
						handle = open(fname, 'w')
						handle.write(r.stdout)
						handle.close()
						assembly = fname
					except Exception as e:
						logger.error("Error during polishing: impossible to write file %s", fname)
						raise e
			else:
				logger.error("Impossible to perform polishing.")

	return assembly


def flye_assemble(args, reads, estimated_genome_size, aslydir):
	'''
	Assemble genome with Flye using long reads
	Returns path to genome assembly or None
	'''
	if bool(estimated_genome_size) is False:
		logger.critical("Impossible to run flye without genome size estimation!")
		sys.exit(1)

	cmd = ["flye", "-o", aslydir, "-g", str(estimated_genome_size), "-t", str(args.threads)]
	if "nanopore" in reads.keys():
		cmd += ["--nano-raw", reads['nanopore']]
	elif "pacbio" in reads.keys():
		cmd += ["--pacbio-raw", reads['pacbio']]

	if args.skip_flye_long_polish:
		cmd += ["--stop-after", "contigger"]

	if run_external(args, cmd).returncode != 0:
		logger.error("Genome assembly finished with errors.")
		logger.critical("Plese check %s for more information.", os.path.join(aslydir, "flye.log"))
		raise Exception("Extermal software error")
	else:
		logger.debug("Assembling finished")
		if args.skip_flye_long_polish:
			assembly = os.path.join(aslydir, "30-contigger/contigs.fasta")
		else:
			assembly = os.path.join(aslydir, "assembly.fasta")
		return assembly


def assemble(args, reads, estimated_genome_size):
	logger.info("Assembling started")
	aslydir = os.path.join(args.output_dir, "assembly")

	if args.assembler == "flye":
		return flye_assemble(args, reads, estimated_genome_size, aslydir)

	if args.assembler == "spades":

		cmd = ["spades.py", "-o", aslydir, "--careful", "-t", str(args.threads),
			"--cov-cutoff", "auto"]

		# Check Spades and get it's version to use verion-specific features
		# as "--merged", "--isolate" etc.
		try:
			spades_stdout = str(subprocess.run(["spades.py", "-v"], stdout=subprocess.PIPE,
				stderr=subprocess.PIPE, universal_newlines=True).stdout)
			spades_version = re.search(r'[\d\.]+', spades_stdout)[0]
			logger.debug("SPAdes version: %s", spades_version)
		except Exception as e:
			logger.critical("Failed to run SPAdes!")
			raise e

		'''
		v_major = int(spades_version.split(".")[0])
		v_minor = int(spades_version.split(".")[1])

		I'm not sure is it worth to use this mode...
		if v_major > 3 or (v_major >= 3 and v_minor >= 14):
			cmd += ['--isolate']
		else:
			cmd += ['--careful']
		'''

		if args.memory_limit:
			cmd += ["-m", str(args.memory_limit)]
		if 'pe' in reads.keys():
			cmd += ["-1", reads['pe'][0], "-2", reads['pe'][1]]
		if 'merged' in reads.keys():
			cmd += ["--merged", reads['merged']]
		if 'single' in reads.keys():
			cmd += ["-s", reads['single']]
		if 'mp' in reads.keys():
			cmd += ["--mp1-1", reads['mp'][0], "--mp1-2", reads['mp'][1]]
		if 'nanopore' in reads.keys():
			cmd += ["--nanopore", reads['nanopore']]
		if 'pacbio' in reads.keys():
			cmd += ["--pacbio", reads['pacbio']]
		if args.no_correction:
			cmd += ["--only-assembler"]
		if args.spades_k_list:
			cmd += ["-k", args.spades_k_list]

		if run_external(args, cmd).returncode != 0:
			logger.error("Genome assembly finished with errors.")
			logger.error("Plese check %s for more information.", os.path.join(aslydir, "spades.log"))
			raise Exception("Extermal software error")
		else:
			logger.debug("Assembling finished")
			if args.use_scaffolds:
				return os.path.join(aslydir, "scaffolds.fasta")
			else:
				return os.path.join(aslydir, "contigs.fasta")

	elif args.assembler == "unicycler":
		# Only to check if able to run unicycler now
		try:
			version_stdout = subprocess.run(["unicycler", "--version"], encoding="utf-8",
				stderr=subprocess.PIPE, stdout=subprocess.PIPE).stdout
			version = re.search(r'v(\d\S*)', version_stdout)[1]
			logger.debug("Unicycler version %s available.", version)
		except Exception as e:
			logger.critical("Failed to run Unicycler!")
			raise e

		cmd = ["unicycler", "-o", aslydir, "-t", str(args.threads),
			"--mode", args.unicycler_mode]
		if 'pe' in reads.keys():
			cmd += ["-1", reads['pe'][0], "-2", reads['pe'][1]]
		if args.no_correction:
			cmd += ["--no_correct"]
		if 'merged' in reads.keys():
			cmd += ["-s", reads['merged']]
		elif 'single' in reads.keys():
			cmd += ["-s", reads['single']]
		if 'nanopore' in reads.keys():
			cmd += ["-l", reads['nanopore']]
		if 'pacbio' in reads.keys():
			cmd += ["-l", reads['pacbio']]

		if run_external(args, cmd).returncode != 0:
			logger.error("Genome assembly finished with errors.")
			logger.error("Plese check %s for more information.", os.path.join(aslydir, "unicycler.log"))
			raise Exception("Extermal software error")
		else:
			logger.debug("Assembling finished")
			assembly = os.path.join(aslydir, "assembly.fasta")
			if args.extract_replicons:
				extract_replicons(args, aslydir)
			return assembly

	else:
		logger.critical("Not yet implemented")
		return None


def extract_replicons(args, aslydir):
	'''
	Unicycler log is not universal for different input files.
	Current implementation works only for shor-read based assemblies.
	'''
	logfile = os.path.join(aslydir, "unicycler.log")
	assemblyfile = os.path.join(aslydir, "assembly.fasta")
	with open(logfile, "r") as log:
		regexp = re.compile(r'\s?\d+.+\scomplete')
		repl_lengths = []
		for line in log.readlines():
			if regexp.search(line):
				_, segments, _, length, _, _, _ = line.split()
				if int(segments) == 1:
					repl_lengths.append(int(length.replace(",", "")))
		if len(repl_lengths) == 0:
			logger.info("No complete replicons found.")
			return 0
		with open(assemblyfile, "r") as assembly:
			replicons = [x for x in SeqIO.parse(assembly, "fasta") if len(x) in repl_lengths]
			for n, replicon in enumerate(replicons, start=1):
				with open(os.path.join(aslydir, f"replicon.{n}.fasta"), "w") as F:
					SeqIO.write(replicon, F, "fasta")
			logger.info("Extracted %s replicon(s).", str(len(repl_lengths)))
			return len(replicons)


def locus_tag_gen(genome):
	logger.info("No locus tag provided. Generating it as MD5 hash of genome")
	with open(genome, 'rb') as genomefile:
		digest = hashlib.md5(genomefile.read()).hexdigest()
		locus_tag = "".join([chr(65 + (int(digest[x], 16) + int(digest[x + 1], 16)) % 26) for x in range(0, 12, 2)])
		logger.info("Locus tag generated: %s", locus_tag)
		return locus_tag


def annotate(args):
	logger.info("Genome annotation started")

	try:
		version_stdout = subprocess.run(["dfast", "--version"], encoding="utf-8",
			stderr=subprocess.PIPE, stdout=subprocess.PIPE).stdout
		version = re.search(r'ver. (\d\S*)', version_stdout)[1]
		logger.debug("DFAST version %s available.", version)
	except Exception as e:
		logger.critical("Failed to run DFAST!")
		raise e

	annodir = os.path.join(args.output_dir, "annotation")

	cmd = ["dfast", "-g", args.genome, "-o", annodir, "--organism",
		" ".join([args.genus, args.species]), "--cpu", str(args.threads)]
	if args.strain:
		cmd += ["--strain", args.strain]
	if args.center_name:
		cmd += ["--center_name", args.center_name]
	if not args.locus_tag:
		args.locus_tag = locus_tag_gen(args.genome)
	if args.locus_tag:
		cmd += ["--locus_tag_prefix", args.locus_tag]
	if args.locus_tag_inc:
		cmd += ["--step", str(args.locus_tag_inc)]
	if args.minimum_contig_length:
		cmd += ["--minimum_length", args.minimum_contig_length]

	if run_external(args, cmd).returncode == 0:
		args.genome = os.path.join(annodir, "genome.fna")
	return args.genome


def check_phix(args):
	logger.info("Checking assembly for presence of Illumina phiX control")
	phix_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data/phiX174.fasta")
	blast_format = "6 sseqid pident slen length"
	cmd = ["blastn", "-query", phix_path, "-subject", args.genome, "-outfmt", blast_format, "-evalue", "1e-6"]

	logger.debug("Running: %s", " ".join(cmd))
	try:
		blast_out = str(subprocess.run(cmd, universal_newlines=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE).stdout).rstrip()
	except Exception as e:
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

	return args.genome


def run_checkm(args):

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

		cmd = ["checkm", "taxonomy_wf", "-f", checkm_outfile, "-x", checkm_ext, "-t", str(args.threads)]

		try:
			checkm_taxon_list = subprocess.run(["checkm", "taxon_list"], encoding="utf-8",
				universal_newlines=True, stderr=subprocess.DEVNULL, stdout=subprocess.PIPE).stdout
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
			args.checkm_taxon = "Archaea" if args.domain == "archaea" else "Bacteria"
			args.checkm_rank = "domain"

		logger.info("%s marker set will be used for CheckM", args.checkm_taxon)

		cmd += [args.checkm_rank, args.checkm_taxon, checkm_indir, checkm_outdir]

	else:
		cmd = ["checkm", "lineage_wf", "-f", checkm_outfile, "-t", str(args.threads), "-x", checkm_ext]
		if not args.checkm_full_tree:
			cmd += ["--reduced_tree"]
		cmd += ["--pplacer_threads", str(args.threads), checkm_indir, checkm_outdir]

	rc = run_external(args, cmd).returncode

	# Cleaning after CheckM
	shutil.rmtree(checkm_indir)
	shutil.rmtree(checkm_outdir)

	if rc == 0:
		return checkm_outfile
	else:
		return None


def check_last_step(args, step):
	if args.last_step == step:
		logger.info("Workflow finished!")
		sys.exit(0)


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

	args.output_dir = os.path.abspath(args.output_dir)
	if os.path.isdir(args.output_dir):
		if args.force:
			shutil.rmtree(args.output_dir)
		else:
			logger.critical("Output directory \"%s\" already exists. Use --force to overwrite.",
				args.output_dir)
			raise FileExistsError(args.output_dir)
	try:
		os.mkdir(args.output_dir)
	except Exception as e:
		logger.critical("Imposible to create directory \"%s\"!", args.output_dir)
		raise e

	fh = logging.FileHandler(os.path.join(args.output_dir, "zga.log"))
	fh.setLevel(logging.DEBUG)
	fh.setFormatter(formatter)
	logger.addHandler(fh)

	logger.info("Pipeline started")

	steps = {"readqc": 1, "processing": 2, "assembling": 3, "polishing": 4,
		"annotation": 6, "check_genome": 5}
	args.first_step = steps[args.first_step]
	args.last_step = steps[args.last_step]
	if args.first_step <= 4:
		reads = check_reads(args)
		logger.debug("Input reads: %s", reads)
		if len(list(reads)) == 0:
			logger.critical("No reads provided for genome assembly")
			raise Exception("No reads provided for genome assembly")
		short_reads_keys = {'pe', 'merged', 'single'}
		short_reads_list = [ reads[k] for k in list(set(reads.keys()) & short_reads_keys) ]

	if args.first_step > 3 and (not args.genome or not os.path.isfile(args.genome)):
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
		if args.genome_size_estimation:
			estimated_genome_size = mash_estimate(args, reads)
		else:
			estimated_genome_size = args.estimated_genome_size
		args.genome = assemble(args, reads, estimated_genome_size)
		check_last_step(args, 3)

	# Short read polishing is only meaningful for flye assembly
	if args.first_step <= 4 and args.perform_polishing:
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
				completeness, contamination, heterogeneity = list(map(float, result.readlines()[3].split()[-3::1]))
				if completeness < 80.0 or (completeness - 5.0 * contamination < 50.0):
					logger.info("The genome assembly has low quality!")
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
