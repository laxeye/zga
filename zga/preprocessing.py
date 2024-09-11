'''This module performs sequencing data preprocessing with BBTools:
filtering, trimming, merging, normalizing.
'''
import logging
import os.path
import shutil
import sys
from zga.essential import create_subdir, run_external


def remove_intermediate(path, *files):
	'''Remove files from directory (path) and keep initial user files'''
	logger = logging.getLogger("main")
	for f in files:
		if os.path.dirname(f) == path and os.path.exists(f):
			logger.debug("Removing %s", f)
			os.remove(f)


def filter_by_tile(args, reads, readdir):
	'''Run filterbytile.sh (BBmap) for Illumina read filtering'''
	logger = logging.getLogger("main")
	for index, lib in enumerate(reads, start=1):
		if (lib["type"] == "short"
			and "forward" in lib.keys()
			and "reverse" in lib.keys()
		):
			filtered_pe_r1 = os.path.join(readdir, f"lib{index}.filtered.r1.fq.gz")
			filtered_pe_r2 = os.path.join(readdir, f"lib{index}.filtered.r2.fq.gz")

			cmd = ["filterbytile.sh", f"in={lib['forward']}", f"in2={lib['reverse']}",
				f"out={filtered_pe_r1}", f"out2={filtered_pe_r2}",
				f"-Xmx={args.memory_limit}G"]

			if run_external(args, cmd):
				remove_intermediate(readdir, lib['forward'], lib['reverse'])
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
			u1 = os.path.join(readdir, f"lib{index}.u1.fq")
			u2 = os.path.join(readdir, f"lib{index}.u2.fq")
			merged = os.path.join(readdir, f"lib{index}.merged.fq")

			cmd = ["bbmerge.sh", "strict=t",
			f"Xmx={args.memory_limit}G", f"t={args.threads}",
			f"in={lib['forward']}", f"in2={lib['reverse']}",
			f"outu1={u1}", f"outu2={u2}", f"out={merged}"]

			if args.bbmerge_extra:
				cmd += args.bbmerge_extra

			if run_external(args, cmd):
				remove_intermediate(readdir, lib['forward'], lib['reverse'])
				lib['forward'], lib['reverse'] = u1, u2
				lib['merged'] = merged

	return reads


def repair_pair(args, readdir, lib, index):
	'''Run repair.sh from BBmap

	Returns:
	(Fixed R1, Fixed R2), Singletons
	'''
	logger = logging.getLogger("main")
	singletons = os.path.join(readdir, f"lib{index}.singletons.fq")
	fixed = (os.path.join(readdir, f"lib{index}.repaired.r1.fq"),
		os.path.join(readdir, f"lib{index}.repaired.r2.fq"))
	cmd = ["repair.sh", f"in={lib[0]}", f"in2={lib[1]}",
		f"out={fixed[0]}", f"out2={fixed[1]}", f"outs={singletons}",
		f"Xmx={args.memory_limit}G"]
	if run_external(args, cmd):
		return (fixed, singletons)
	logger.error("Error during repair of paired-end reads %s and %s",
		lib[0], lib[1])
	sys.exit(1)


def bbduk_process(args, reads, readdir):
	"""Perform trimming and filtering of short reads"""
	logger = logging.getLogger("main")
	precmd = ["bbduk.sh", f"Xmx={args.memory_limit}G", f"t={args.threads}",
			f"ref={args.adapters}", f"k={args.bbduk_k}", "ktrim=r",
			"qtrim=rl", f"trimq={args.quality_cutoff}",
			f"entropy={args.entropy_cutoff}",
			f"minlength={args.min_short_read_length}"]
	if args.bbduk_extra:
		precmd += args.bbduk_extra

	for index, lib in enumerate(
		[lib for lib in reads if lib["type"] == "short"],
		start=1):

		if "forward" in lib.keys() and "reverse" in lib.keys():
			logger.info("Trimming and filtering paired end reads")
			initial = (lib['forward'], lib['reverse'])
			out_pe1 = os.path.join(readdir, f"lib{index}.r1.fq")
			out_pe2 = os.path.join(readdir, f"lib{index}.r2.fq")
			out_stats = os.path.join(readdir, f"lib{index}.bbduk.pe.txt")
			cmd = precmd + [f"in={lib['forward']}", f"in2={lib['reverse']}",
				f"out={out_pe1}", f"out2={out_pe2}", f"stats={out_stats}"]

			if run_external(args, cmd):
				remove_intermediate(readdir, *initial)
				lib['forward'], lib['reverse'] = out_pe1, out_pe2
			else:
				logger.error(
					"Error during processing paired-end reads %s and %s",
					lib['forward'], lib['reverse']
				)
				logger.warning("Trying to repair paired-end reads.")
				fixed, discarded = repair_pair(args, readdir, initial, index)
				remove_intermediate(readdir, *initial)
				cmd = precmd + [f"in={fixed[0]}", f"in2={fixed[1]}",
					f"out={out_pe1}", f"out2={out_pe2}", f"stats={out_stats}"]
				if run_external(args, cmd):
					remove_intermediate(readdir, *fixed)
					lib['forward'], lib['reverse'] = out_pe1, out_pe2
					if 'single' not in lib.keys():
						lib['single'] = discarded

		for read_type in ["single", "merged"]:
			if read_type in lib.keys():
				logger.info("Trimming and filtering %s reads", read_type)
				out = os.path.join(readdir, f"lib{index}.{read_type}.fq")
				out_stats = os.path.join(readdir, f"lib{index}.bbduk.{read_type}.txt")
				cmd = precmd + [f"in={lib[read_type]}", f"out={out}",
					f"stats={out_stats}"]

				if run_external(args, cmd):
					remove_intermediate(readdir, lib[read_type])
					lib[read_type] = out

	return reads


def tadpole_correct(args, reads, readdir):
	'''Correct short reads with tadpole'''
	logger = logging.getLogger("main")
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

			if run_external(args, cmd):
				remove_intermediate(readdir, *initial)
				lib['forward'], lib['reverse'] = out_pe1, out_pe2

		for read_type in ["single", "merged"]:
			if read_type in lib.keys():
				logger.info("Error correction of %s reads", read_type)
				initial = lib[read_type]
				out = os.path.join(readdir, f"lib{index}.ecc.{read_type}.fq")
				cmd = precmd + [f"in={initial}", f"out={out}"]

				if run_external(args, cmd):
					remove_intermediate(readdir, initial)
					lib[read_type] = out

	return reads


def bbnorm(args, reads, readdir):
	'''Normalize reads by k-mer coverage depth with BBnorm'''
	logger = logging.getLogger("main")
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

			if run_external(args, cmd):
				remove_intermediate(readdir, *initial)
				lib['forward'], lib['reverse'] = out_pe1, out_pe2

		for read_type in ["single", "merged"]:
			if read_type in lib.keys():
				logger.info("Normalization of %s reads", read_type)
				initial = lib[read_type]
				out = os.path.join(readdir, f"lib{index}.norm.{read_type}.fq")
				cmd = precmd + [f"in={initial}", f"out={out}"]

				if run_external(args, cmd):
					remove_intermediate(readdir, initial)
					lib[read_type] = out

	return reads


def compress_reads(args, reads, readdir):
	'''Compress reads with pigz or gzip after processing'''
	logger = logging.getLogger("main")
	if shutil.which('pigz'):
		cmd = ['pigz', '-p', str(args.threads)]
	else:
		cmd = ['gzip']

	for lib in reads:
		for k, v in lib.items():
			f = os.path.join(readdir, v)
			if os.path.isfile(f) and os.path.splitext(f)[-1] not in ('gz', 'bz2'):
				logger.debug("Compressing %s", v)
				if run_external(args, cmd + [f]):
					lib[k] = f'{v}.gz'

	return reads


def mp_read_processing(args, reads, readdir):
	'''Filtering Illumina mate pair reads'''
	logger = logging.getLogger("main")
	for index, lib in enumerate([lib for lib in reads if lib["type"] == "mate-pair"], start=1):
		prefix = os.path.join(readdir, f"mate.{index}")

		cmd = [
			"nxtrim", "-1", lib['forward'], "-2", lib['reverse'],
			"--separate", "--rf", "--justmp", "-O", prefix, "-l",
			str(args.min_short_read_length)
		]
		logger.info("Processing mate-pair reads.")
		if run_external(args, cmd):
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


def read_processing(args, reads):
	'''Pipeline for read processing

	Returns
	reads (dict)
	'''
	logger = logging.getLogger("main")
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
	if args.bbmerge or args.bbmerge_extra:
		reads = merge_bb(args, reads, readdir)

	# Processing Illumina mate-pairs
	if not args.no_nxtrim:
		reads = mp_read_processing(args, reads, readdir)

	reads = compress_reads(args, reads, readdir)

	logger.info("Read processing finished")

	return reads
