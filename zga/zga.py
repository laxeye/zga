#!/usr/bin/env python3
import argparse
import os.path
import logging
import sys
import shutil
import subprocess
import re
import hashlib
import json
from Bio import SeqIO
from zga import __version__
from zga.assemblers import assemble


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
	general_args.add_argument("-t", "--threads", type=int, default=1,
		help="Number of CPU threads to use (where possible)")
	general_args.add_argument("-m", "--memory-limit", type=int, default=8,
		help="Memory limit in GB (default 8)")
	general_args.add_argument("--genus", default="Unknown",
		help="Provide genus if known")
	general_args.add_argument("--species", default="sp.",
		help="Provide species if known")
	general_args.add_argument("--strain",
		help="Provide strain if known")
	general_args.add_argument("--transparent", action="store_true",
		help="Show output from tools inside the pipeline")
	general_args.add_argument("--domain",
		default="bacteria", choices=['archaea', 'bacteria'],
		help="Provide prokaryotic domain: bacteria or archaea")
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
	reads_args.add_argument("--bbduk-extra", nargs='*',
		help="Extra options for BBduk. Should be space-separated.")
	reads_args.add_argument("--tadpole-correct", action="store_true",
		help="Perform error correction of short reads with tadpole.sh from BBtools."
		+ "SPAdes correction may be disabled with \"--no-spades-correction\".")
	reads_args.add_argument("--bbmerge", action="store_true",
		help="Merge overlapped paired-end reads with BBMerge.")
	reads_args.add_argument("--bbmerge-extend", type=int,
		help="Perform k-mer read extension by specified length "
		+ "if initial merging wasn't succesfull.")
	reads_args.add_argument("--bbmerge-extend-kmer", type=int, default=40,
		help="K-mer length for read extension, default 40.")
	reads_args.add_argument("--bbmerge-trim", type=int,
		help="Before merging trim bases with phred score less than a specified value.")
	reads_args.add_argument("--normalize-kmer-cov", type=int,
		help="Normalize read depth based on kmer counts to arbitrary value.")
	reads_args.add_argument("--calculate-genome-size", action="store_true",
		help="Estimate genome size with mash.")
	reads_args.add_argument("--genome-size-estimation", type=int,
		help="Genome size in bp (no K/M suffix supported) for Flye assembler, if known.")
	reads_args.add_argument("--mash-kmer-copies", type=int, default=10,
		help="Minimum copies of each k-mer to include in size estimation")
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
		help="Number of polishing iterations.")

	check_args = parser.add_argument_group(title="Genome check settings")
	# phiX
	check_args.add_argument("--check-phix", action="store_true",
		help="Check genome for presence of PhiX control sequence.")
	# CheckM
	check_args.add_argument("--checkm-mode",
		default="taxonomy_wf", choices=['taxonomy_wf', 'lineage_wf'],
		help="Select CheckM working mode. Default is checking for domain-specific marker-set.")
	check_args.add_argument("--checkm-rank",
		help="Rank of taxon for CheckM. Run 'checkm taxon_list' for details.")
	check_args.add_argument("--checkm-taxon",
		help="Taxon for CheckM. Run 'checkm taxon_list' for details.")
	check_args.add_argument("--checkm-full-tree", action="store_true",
		help="Use full tree for inference of marker set, requires LOTS of memory.")

	anno_args = parser.add_argument_group(title="Annotation settings")
	anno_args.add_argument("-g", "--genome",
		help="Genome assembly (when starting from annotation).")
	anno_args.add_argument("--gcode", default=11, type=int,
		help="Genetic code.")
	anno_args.add_argument("--locus-tag",
		help="Locus tag prefix. If not provided prefix will be generated from MD5 checksum.")
	anno_args.add_argument("--locus-tag-inc", default=10, type=int,
		help="Locus tag increment, default = 10")
	anno_args.add_argument("--center-name", help="Genome center name.")
	anno_args.add_argument("--minimum-contig-length",
		help="Minimum sequence length in genome assembly.")
	anno_args.add_argument("--dfast-config",
		help="Custom DFAST configuration file.")

	args = parser.parse_args()

	if (args.assembler == 'spades'
		and not (
			bool(args.pe_1)
			or bool(args.pe_merged)
			or bool(args.single_end)
			or bool(args.mp_1)
		)
	):
		logger.error("Impossible to run SPAdes without short reads!")
		raise Exception("Bad parameters.")

	if args.assembler == 'flye':

		if args.nanopore and args.pacbio:
			logger.error("Impossible to run Flye on mixed long reads!")
			raise Exception("Bad parameters.")

		if not (bool(args.nanopore) or bool(args.pacbio)):
			logger.error("Impossible to run Flye without long reads!")
			raise Exception("Bad parameters.")

		if not args.genome_size_estimation:
			args.calculate_genome_size = True
			logger.info("Genome size was not provided. It will be calculated with mash.")

	if (args.flye_short_polish
		or args.first_step == 'polishing'
		or args.last_step == 'polishing'
	):
		args.perform_polishing = True

	if args.dfast_config:
		if os.path.isfile(args.dfast_config):
			args.dfast_config = os.path.abspath(args.dfast_config)
		else:
			logger.error("File \"%s\" not found!", args.dfast_config)
			raise FileNotFoundError(
				"DFAST config file \"%s\" not found." % args.dfast_config
			)

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
		for f in reads:
			if f != "n/a" and not os.path.isfile(f):
				logger.error("File %s doesn't exist", f)
				raise FileNotFoundError("File %s doesn't exist" % f)

	short_libs = {}
	if args.pe_1 and args.pe_2:
		short_libs["forward"] = args.pe_1
		short_libs["reverse"] = args.pe_2
	if args.single_end:
		short_libs["single"] = args.single_end
	if args.pe_merged:
		short_libs["merged"] = args.pe_merged

	if len(short_libs.values()) > 0:
		short_lib_N = max(list(map(len, short_libs.values())))
		for i in range(short_lib_N):
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


def create_subdir(parent, child) -> str:
	'''Tries to create a subdirectory, returns path if succes.'''
	path = os.path.join(parent, child)
	try:
		os.mkdir(path)
	except Exception as e:
		logger.critical("Impossible to create directory \"%s\"", path)
		raise e
	return path


def run_external(args, cmd, keep_stdout=False, keep_stderr=True):
	'''Run external command using subprocess.

	Returns subprocess.CompletedProcess
	'''
	logger.debug("Running: %s", " ".join(cmd))
	stderr_dest = subprocess.PIPE if keep_stderr else subprocess.DEVNULL
	stdout_dest = subprocess.PIPE if keep_stdout else subprocess.DEVNULL

	try:
		r = subprocess.run(cmd, check=True, stderr=stderr_dest, stdout=stdout_dest, encoding="utf-8")
		if args.transparent:
			print(r.stderr, file=sys.stderr)
		return r
	except subprocess.CalledProcessError as e:
		logger.error("Error during execution of: %s", ' '.join(e.cmd))
		logger.info("Please see the logfile for additional information: %s",
			os.path.join(args.output_dir, "zga.log"))
		logger.debug("External tool stderr:\n%s", e.stderr)
		return None


def read_qc(args, reads):
	'''Perform read QC with fastp'''
	logger.info("Read quality control started")
	qcoutdir = create_subdir(args.output_dir, "readQC")
	precmd = ["fastp", "-L", "-Q", "-G", "-A", "-z", "1", "--stdout",
		"-w", str(args.threads)]
	for lib in reads:
		for t, r in lib.items():
			if t == "type":
				continue
			prefix = os.path.join(qcoutdir, os.path.split(r)[-1])
			cmd = precmd + ["-i", r, "-h", f"{prefix}.html", "-j", f"{prefix}.json"]
			logger.debug("QC of %s", r)
			run_external(args, cmd)


def remove_intermediate(path, *files):
	'''Remove files from directory (path) and keep initial user files'''
	for f in files:
		if os.path.dirname(f) == path and os.path.exists(f):
			logger.debug("Removing %s", f)
			os.remove(f)


def filter_by_tile(args, reads, readdir):
	'''Run filterbytile.sh (BBmap) for Illumina read filtering'''
	for index, lib in enumerate(reads, start=1):
		if (lib["type"] == "short"
			and "forward" in lib.keys()
			and "reverse" in lib.keys()
		):
			initial = (lib['forward'], (lib['reverse']))
			filtered_pe_r1 = os.path.join(readdir, f"lib{index}.filtered.r1.fq.gz")
			filtered_pe_r2 = os.path.join(readdir, f"lib{index}.filtered.r2.fq.gz")

			cmd = ["filterbytile.sh", f"in={initial[0]}", f"in2={initial[1]}",
				f"out={filtered_pe_r1}", f"out2={filtered_pe_r2}",
				f"-Xmx={args.memory_limit}G"]

			if run_external(args, cmd) is not None:
				remove_intermediate(readdir, *initial)
				lib['forward'], lib['reverse'] = filtered_pe_r1, filtered_pe_r2
			else:
				logger.warning("Filtering by tile wasn't perfomed correctly.")

	return reads


def merge_bb(args, reads, readdir):
	'''Performs merging (and extension) of overlapping paired-end reads.

	Parameters:
	args (argparse.Namespace)
	reads (dict) : input reads
	readdir (path) : path to reads output directory

	Returns:
	reads (dict) : modified if bbmerge returns 0.
	'''
	for index, lib in enumerate(reads, start=1):
		if (lib["type"] == "short"
			and "forward" in lib.keys()
			and "reverse" in lib.keys()
			and "merged" not in lib.keys()
		):
			# Initial, unmerged and merged filenames
			initial = (lib['forward'], (lib['reverse']))
			u1 = os.path.join(readdir, f"lib{index}.u1.fq")
			u2 = os.path.join(readdir, f"lib{index}.u2.fq")
			merged = os.path.join(readdir, f"lib{index}.merged.fq")

			cmd = ["bbmerge.sh", f"Xmx={args.memory_limit}G", f"t={args.threads}",
			f"in={initial[0]}", f"in2={initial[1]}",
			f"outu1={u1}", f"outu2={u2}", f"out={merged}"]

			if args.bbmerge_extend or bbmerge_extend_kmer:
				cmd += [f"extend2={args.bbmerge_extend}",
				f"k={args.bbmerge_extend_kmer}", "rsem=t"]
			else:
				cmd.append("strict=t")

			if args.bbmerge_trim:
				cmd += ["qtrim2=t", f"trimq={args.bbmerge_trim}"]
			logger.info("Merging paired-end reads.")

			if run_external(args, cmd) is not None:
				remove_intermediate(readdir, *initial)
				lib['forward'], lib['reverse'] = u1, u2
				lib['merged'] = merged

	return reads


def repair_pair(args, readdir, lib, index):
	'''Run repair.sh from BBmap

	Returns:
	(Fixed R1, Fixed R2), Singletons
	'''
	singletons = os.path.join(readdir, f"lib{index}.singletons.fq")
	fixed = (os.path.join(readdir, f"lib{index}.repaired.r1.fq"),
		os.path.join(readdir, f"lib{index}.repaired.r2.fq"))
	cmd = ["repair.sh", f"in={lib[0]}", f"in2={lib[1]}",
		f"out={fixed[0]}", f"out2={fixed[1]}", f"outs={singletons}",
		f"Xmx={args.memory_limit}G"]
	if run_external(args, cmd) is not None:
		return (fixed, singletons)
	logger.error("Error during repair of paired-end reads %s and %s",
		lib[0], lib[1])
	sys.exit(1)


def bbduk_process(args, reads, readdir):
	"""Perform trimming and filtering of short reads"""
	bbduk_kmer = 19  # K-mer length for contaminant/adapter removal
	precmd = ["bbduk.sh", f"Xmx={args.memory_limit}G", f"t={args.threads}",
			f"ref={args.adapters}", f"k={bbduk_kmer}", "ktrim=r",
			"qtrim=r", f"trimq={args.quality_cutoff}",
			f"entropy={args.entropy_cutoff}",
			f"minlength={args.min_short_read_length}"]
	if args.bbduk_extra:
		precmd += bbduk_extra

	for index, lib in enumerate(
		[lib for lib in reads if lib["type"] == "short"],
		start=1):

		if "forward" in lib.keys() and "reverse" in lib.keys():
			logger.info("Trimming and filtering paired end reads")
			initial = (lib['forward'], (lib['reverse']))
			out_pe1 = os.path.join(readdir, f"lib{index}.r1.fq")
			out_pe2 = os.path.join(readdir, f"lib{index}.r2.fq")
			out_stats = os.path.join(readdir, f"lib{index}.bbduk.pe.txt")
			cmd = precmd + [f"in={initial[0]}", f"in2={initial[1]}",
				f"out={out_pe1}", f"out2={out_pe2}", f"stats={out_stats}"]

			if run_external(args, cmd) is not None:
				remove_intermediate(readdir, *initial)
				lib['forward'], lib['reverse'] = out_pe1, out_pe2
			else:
				logger.error(
					"Error during processing paired-end reads %s and %s",
					initial[0], initial[1]
				)
				logger.warning("Trying to repair paired-end reads.")
				fixed, discarded = repair_pair(args, readdir, initial, index)
				remove_intermediate(readdir, *initial)
				cmd = precmd + [f"in={fixed[0]}", f"in2={fixed[1]}",
					f"out={out_pe1}", f"out2={out_pe2}", f"stats={out_stats}"]
				if run_external(args, cmd) is not None:
					remove_intermediate(readdir, *fixed)
					lib['forward'], lib['reverse'] = out_pe1, out_pe2
					if 'single' not in lib.keys():
						lib['single'] = discarded

		for read_type in ["single", "merged"]:
			if read_type in lib.keys():
				logger.info("Trimming and filtering %s reads", read_type)
				initial = lib[read_type]
				out = os.path.join(readdir, f"lib{index}.{read_type}.fq")
				out_stats = os.path.join(readdir, f"lib{index}.bbduk.{read_type}.txt")
				cmd = precmd + [f"in={initial}", f"out={out}",
					f"stats={out_stats}"]

				if run_external(args, cmd) is not None:
					remove_intermediate(readdir, initial)
					lib[read_type] = out

	return reads


def tadpole_correct(args, reads, readdir):
	'''Correct short reads with tadpole'''
	precmd = ["tadpole.sh",
		f"Xmx={args.memory_limit}G",
		f"t={args.threads}",
		"mode=correct"]

	for index, lib in enumerate([lib for lib in reads if lib["type"] == "short"], start=1):
		if "forward" in lib.keys() and "reverse" in lib.keys():
			logger.info("Error correction of paired end reads")
			initial = (lib['forward'], (lib['reverse']))
			out_pe1 = os.path.join(readdir, f"lib{index}.ecc.pe.r1.fq")
			out_pe2 = os.path.join(readdir, f"lib{index}.ecc.pe.r2.fq")
			cmd = precmd + [f"in={initial[0]}", f"in2={initial[1]}",
				f"out={out_pe1}", f"out2={out_pe2}"]

			if run_external(args, cmd) is not None:
				remove_intermediate(readdir, *initial)
				lib['forward'], lib['reverse'] = out_pe1, out_pe2

		for read_type in ["single", "merged"]:
			if read_type in lib.keys():
				logger.info("Error correction of %s reads", read_type)
				initial = lib[read_type]
				out = os.path.join(readdir, f"lib{index}.ecc.{read_type}.fq")
				cmd = precmd + [f"in={initial}", f"out={out}"]

				if run_external(args, cmd) is not None:
					remove_intermediate(readdir, initial)
					lib[read_type] = out

	return reads


def bbnorm(args, reads, readdir):
	'''Normalize reads by k-mer coverage depth with BBnorm'''
	precmd = ["bbnorm.sh",
		f"Xmx={args.memory_limit}G",
		f"t={args.threads}",
		f"minq={args.quality_cutoff}",
		f"target={args.normalize_kmer_cov}"]

	for index, lib in enumerate([lib for lib in reads if lib["type"] == "short"], start=1):
		if "forward" in lib.keys() and "reverse" in lib.keys():
			logger.info("Normalization of paired end reads")
			initial = (lib['forward'], (lib['reverse']))
			out_pe1 = os.path.join(readdir, f"lib{index}.norm.pe.r1.fq")
			out_pe2 = os.path.join(readdir, f"lib{index}.norm.pe.r2.fq")
			cmd = precmd + [f"in={initial[0]}", f"in2={initial[1]}",
				f"out={out_pe1}", f"out2={out_pe2}"]

			if run_external(args, cmd) is not None:
				remove_intermediate(readdir, *initial)
				lib['forward'], lib['reverse'] = out_pe1, out_pe2

		for read_type in ["single", "merged"]:
			if read_type in lib.keys():
				logger.info("Normalization of %s reads", read_type)
				initial = lib[read_type]
				out = os.path.join(readdir, f"lib{index}.norm.{read_type}.fq")
				cmd = precmd + [f"in={initial}", f"out={out}"]

				if run_external(args, cmd) is not None:
					remove_intermediate(readdir, initial)
					lib[read_type] = out

	return reads


def compress_reads(args, reads, readdir):
	'''Compress reads with pigz or gzip after processing'''

	if shutil.which('pigz') is not None:
		cmd = ['pigz', '-p', str(args.threads)]
	else:
		cmd = ['gzip']

	for lib in reads:
		for k, v in lib.items():
			f = os.path.join(readdir, v)
			if os.path.isfile(f) and os.path.splitext(f)[-1] not in ('gz', 'bz2'):
				logger.debug("Compressing %s", v)
				if run_external(args, cmd + [f]) is not None:
					lib[k] = f'{v}.gz'

	return reads


def read_processing(args, reads):
	'''Pipeline for read processing

	Returns
	reads (dict)
	'''
	logger.info("Reads processing started")
	readdir = create_subdir(args.output_dir, "reads")
	sr_adapters = os.path.join(
		os.path.dirname(os.path.abspath(__file__)),
		"data/sr.adapters.fasta")

	if args.adapters and os.path.isfile(args.adapters):
		args.adapters = os.path.abspath(args.adapters)
	else:
		args.adapters = sr_adapters

	if args.filter_by_tile:
		reads = filter_by_tile(args, reads, readdir)

	reads = bbduk_process(args, reads, readdir)

	if args.tadpole_correct:
		reads = tadpole_correct(args, reads, readdir)

	if args.normalize_kmer_cov:
		reads = bbnorm(args, reads, readdir)

	# Merging overlapping paired-end reads
	if args.bbmerge or args.bbmerge_trim or bbmerge_extend:
		reads = merge_bb(args, reads, readdir)

	# Processing Illumina mate-pairs
	if not args.no_nxtrim:
		reads = mp_read_processing(args, reads, readdir)

	reads = compress_reads(args, reads, readdir)

	logger.info("Read processing finished")

	return reads


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

	if r is not None:
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
	'''Filtering Illumina mate pair reads'''
	for index, lib in enumerate([lib for lib in reads if lib["type"] == "mate-pair"], start=1):
		prefix = os.path.join(readdir, f"mate.{index}")

		cmd = [
			"nxtrim", "-1", lib['forward'], "-2", lib['reverse'],
			"--separate", "--rf", "--justmp", "-O", prefix, "-l",
			str(args.min_short_read_length)
		]
		logger.info("Processing mate-pair reads.")
		if run_external(args, cmd) is not None:
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
				lib['forward'], lib['reverse'] = (
					f"{prefix}_R1.all.fastq.gz", f"{prefix}_R2.all.fastq.gz")
			else:
				lib['forward'], lib['reverse'] = (
					f"{prefix}_R1.mp.fastq.gz", f"{prefix}_R2.mp.fastq.gz")

	return reads


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
	if run_external(args, cmd) is not None:
		return target
	else:
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
				if r is not None:
					digest = hashlib.md5(r.stdout.encode('utf-8')).hexdigest()
					suffix = "".join(
						[chr(65 + (int(digest[x], 16) + int(digest[x + 1], 16)) % 26) for x in range(0, 20, 2)]
					)
					fname = os.path.join(polish_dir, f"polished.{suffix}.fna")
					try:
						handle = open(fname, 'w')
						handle.write(r.stdout)
						handle.close()
						assembly = fname
					except Exception:
						logger.error("Error during polishing: impossible to write file %s", fname)
			else:
				logger.error("Impossible to perform polishing.")

	return assembly


def locus_tag_gen(genome) -> str:
	'''Returns a six letter locus tag generated from MD5 of genome assembly'''
	logger.info("No locus tag provided. Generating it as MD5 hash of genome")
	with open(genome, 'rb') as genomefile:
		digest = hashlib.md5(genomefile.read()).hexdigest()
		locus_tag = "".join(
			[chr(65 + (int(digest[x], 16) + int(digest[x + 1], 16)) % 26) for x in range(0, 12, 2)]
		)
		logger.info("Locus tag generated: %s", locus_tag)
		return locus_tag


def annotate(args) -> str:
	'''Returns path to annotated genome (FASTA)'''
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
	if args.dfast_config:
		cmd += ["--config", args.dfast_config]

	if run_external(args, cmd):
		args.genome = os.path.join(annodir, "genome.fna")
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
			cmd, universal_newlines=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE
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
				stdout=subprocess.PIPE
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
			args.checkm_taxon = "Archaea" if args.domain == "archaea" else "Bacteria"
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

	if r is not None:
		return checkm_outfile
	else:
		logger.error("CheckM didn't finish properly.")
		return None


def get_N_L_metric(lengths, value=50):
	'''Returns a tuple containing NX and LX metric'''
	l_total = sum(lengths)
	metric = 0.01 * value
	l_sum = 0
	for i, x in enumerate(lengths, 1):
		l_sum += x
		if l_sum >= l_total * metric:
			return x, i


def assembly_stats(genome):
	'''Returns a dict containing genome stats'''
	seq_records = SeqIO.parse(genome, "fasta")
	lengths = sorted([len(x.seq) for x in seq_records], reverse=True)
	stats = {
		'Sequence count': len(lengths),
		'Total length': sum(lengths),
		'Max length': lengths[0]
	}
	stats['N50'], stats['L50'] = get_N_L_metric(lengths, 50)
	stats['N90'], stats['L90'] = get_N_L_metric(lengths, 90)
	return stats


def write_assembly_stats(args, stats, prefix, s_format="human"):
	'''Write assembly stats to a file'''
	ext_dict = {"human": "txt", "json": "json", "table": "tsv"}
	filename = f"{prefix}.assembly.{ext_dict[s_format]}"
	with open(os.path.join(args.output_dir, filename), 'w') as dest:
		if s_format == "human":
			for k, v in stats.items():
				print(f"{k}\t{v}", file=dest)
		if s_format == "json":
			print(json.dumps(stats), file=dest)
		if s_format == "table":
			header = "\t".join(stats.keys())
			data = "\t".join(list(map(str, stats.values())))
			print(f"#{header}", file=dest)
			print(f"{data}", file=dest)
	return filename


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
		stats = assembly_stats(args.genome)
		logger.info("Assembly length: %s", stats['Total length'])
		logger.info("Contig count: %s", stats['Sequence count'])
		logger.info("N50: %s", stats['N50'])
		write_assembly_stats(args, stats, prefix=args.assembler, s_format="table")

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
