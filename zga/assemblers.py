'''This module runs genome assembling with one of four assemblers:
Flye, Unicycler, SPAdes or MEGAHIT
'''
import json
import logging
import os.path
import re
import subprocess
import sys

from Bio import SeqIO
from zga.essential import run_external


def get_assembler_version(assembler):
	'''Get version of spades, unicycler or MEGAHIT'''

	logger = logging.getLogger("main")

	try:
		version_stdout = subprocess.run(
			[assembler, "--version"], encoding="utf-8",
			stderr=subprocess.PIPE, stdout=subprocess.PIPE, check=True
		).stdout
		version = re.search(r'v(\d\S*)', version_stdout)[1]
		logger.debug("%s version %s available.", assembler, version)
		return version
	except Exception as ex:
		logger.critical("Failed to run %s!", assembler)
		raise ex


def flye_assemble(args, reads, estimated_genome_size, aslydir) -> str:
	'''Assemble genome with Flye using long reads

	Returns:
	assembly (str) : path to genome assembly
	'''

	logger = logging.getLogger("main")

	if bool(estimated_genome_size) is False:
		logger.critical("Impossible to run flye without genome size estimation!")
		sys.exit(1)

	cmd = ["flye", "-o", aslydir,
		"-g", str(estimated_genome_size),
		"-t", str(args.threads)]

	types = [x["type"] for x in reads]
	if "nanopore" in types:
		long_reads = "nanopore"
		flye_key = "--nano-raw"
	elif "pacbio" in types:
		long_reads = "pacbio"
		flye_key = "--pacbio-raw"
	else:
		raise Exception("No long reads provided for Flye")

	reads = [x['single'] for x in reads if x['type'] == long_reads]
	cmd += [flye_key, *reads]

	if args.flye_skip_long_polish:
		cmd += ["--stop-after", "contigger"]

	if run_external(args, cmd) is not None:
		logger.debug("Assembling finished")
		if args.flye_skip_long_polish:
			assembly = os.path.join(aslydir, "30-contigger/contigs.fasta")
		else:
			assembly = os.path.join(aslydir, "assembly.fasta")
		return assembly

	return assembler_failed(aslydir, "flye.log")


def spades_assemble(args, reads, aslydir) -> str:
	'''Performs genome assembly with SPAdes.

	Returns:
	assembly (str) : path to genome assembly
	'''

	logger = logging.getLogger("main")
	version = get_assembler_version('spades.py')
	cmd = ["spades.py", "-o", aslydir, "--careful", "-t", str(args.threads),
		"--cov-cutoff", "auto", "-m", str(args.memory_limit)]

	for index, lib in enumerate(reads, start=1):
		if lib['type'] == "short":
			if "forward" in lib.keys() and "reverse" in lib.keys():
				cmd += [f"--pe{index}-1", lib['forward'],
					f"--pe{index}-2", lib['reverse']]
			if "merged" in lib.keys():
				cmd += [f"--pe{index}-m", lib['merged']]
			if "single" in lib.keys():
				cmd += [f"--pe{index}-s", lib['single']]
		if lib['type'] == 'mate-pair':
			cmd += [f"--mp{index}-1", lib['forward'],
				f"--mp{index}-2", lib['reverse']]
		if lib['type'] == 'nanopore':
			cmd += ["--nanopore", lib['single']]
		if lib['type'] == 'pacbio':
			cmd += ["--pacbio", lib['single']]

	if args.no_spades_correction:
		cmd += ["--only-assembler"]
	if args.spades_k_list:
		cmd += ["-k", args.spades_k_list]

	if run_external(args, cmd) is not None:
		logger.debug("Assembling with SPAdes %s finished", version)
		if args.use_scaffolds:
			return os.path.join(aslydir, "scaffolds.fasta")
		return os.path.join(aslydir, "contigs.fasta")

	return assembler_failed(aslydir, "spades.log")


def unicycler_assemble(args, reads, aslydir) -> str:
	'''Performs genome assembly with Unicycler.

	Returns:
	assembly (str) : path to assembled genome
	'''

	logger = logging.getLogger("main")
	version = get_assembler_version('unicycler')

	cmd = ["unicycler", "-o", aslydir, "-t", str(args.threads),
		"--mode", args.unicycler_mode,
		"--spades_options", f'--memory {args.memory_limit}']

	for lib in reads:
		sr_parsed = False
		if lib['type'] == "short":
			if sr_parsed:
				continue
			if "forward" in lib.keys() and "reverse" in lib.keys():
				cmd += ["-1", lib['forward'], "-2", lib['reverse']]
			if "merged" in lib.keys():
				cmd += ["-s", lib['merged']]
			elif "single" in lib.keys():
				cmd += ["-s", lib['single']]
			sr_parsed = True
		if lib['type'] == 'nanopore':
			cmd += ["-l", lib['single']]
		if lib['type'] == 'pacbio':
			cmd += ["-l", lib['single']]

	if run_external(args, cmd) is not None:
		logger.info("Assembling with Unicycler %s finished", version)
		assembly = os.path.join(aslydir, "assembly.fasta")
		if args.extract_replicons:
			extract_replicons(aslydir)
		return assembly

	return assembler_failed(aslydir, "unicycler.log")


def extract_replicons(aslydir) -> int:
	'''Extracts complete replicons from short read assemblies.

	Unicycler log is not universal for different input files.
	Current implementation works only for short read based assemblies.

	Returns number of extracted replicons.
	'''
	logger = logging.getLogger("main")
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
			contigs = SeqIO.parse(assembly, "fasta")
			replicons = [x for x in contigs if len(x) in repl_lengths]
			for i, replicon in enumerate(replicons, start=1):
				with open(os.path.join(aslydir, f"replicon.{i}.fasta"), "w") as repl_file:
					SeqIO.write(replicon, repl_file, "fasta")
			logger.info("Extracted %s replicon(s).", str(len(repl_lengths)))
			return len(replicons)


def megahit_assemble(args, reads, aslydir) -> str:
	'''Performs genome assembly with MEGAHIT.

	Returns:
	assembly (str) : path to assembled genome
	'''

	logger = logging.getLogger("main")
	version = get_assembler_version('megahit')

	cmd = ["megahit", "-o", aslydir, "-t", str(args.threads),
		"--memory", str(args.memory_limit * 10**9)]

	se_reads = []
	forward_pe_reads = []
	reverse_pe_reads = []

	for lib in reads:
		if lib['type'] != 'short':
			continue
		if 'single' in lib.keys():
			se_reads.append(lib['single'])
		if 'merged' in lib.keys():
			se_reads.append(lib['merged'])
		if 'forward' in lib.keys() and 'reverse' in lib.keys():
			forward_pe_reads.append(lib['forward'])
			reverse_pe_reads.append(lib['reverse'])

	if len(forward_pe_reads) > 0:
		cmd += ['-1', ','.join(forward_pe_reads)]
		cmd += ['-2', ','.join(reverse_pe_reads)]
	if len(se_reads) > 0:
		cmd += ['-r', ','.join(se_reads)]

	if run_external(args, cmd) is not None:
		logger.debug("Assembling with MEGAHIT %s finished", version)
		return os.path.join(aslydir, "final.contigs.fa")

	return assembler_failed(aslydir, "log")


def assembler_failed(aslydir, logfile):
	'''Logs assembler's fail and raises an exception'''
	logger = logging.getLogger("main")
	logger.error("Genome assembly finished with errors.")
	logger.error("Plese check %s for more information.",
		os.path.join(aslydir, logfile))
	raise Exception("External software error")


def get_nx_lx_metric(lengths, value=50):
	'''Returns a tuple containing NX and LX metric'''
	l_total = sum(lengths)
	metric = 0.01 * value
	l_sum = 0
	for i, len_i in enumerate(lengths, 1):
		l_sum += len_i
		if l_sum >= l_total * metric:
			return i, len_i
	raise Exception("Error during assembly statistics calculation.")


def assembly_stats(genome):
	'''Returns a dict containing genome stats'''
	seq_records = SeqIO.parse(genome, "fasta")
	lengths = sorted([len(x.seq) for x in seq_records], reverse=True)
	stats = {
		'Sequence count': len(lengths),
		'Total length': sum(lengths),
		'Max length': lengths[0]
	}
	stats['L50'], stats['N50'] = get_nx_lx_metric(lengths, 50)
	stats['L90'], stats['N90'] = get_nx_lx_metric(lengths, 90)
	return stats


def write_assembly_stats(args, stats, prefix, s_format="human"):
	'''Write assembly stats to a file'''
	ext_dict = {"human": "txt", "json": "json", "table": "tsv"}
	filename = f"{prefix}.assembly.{ext_dict[s_format]}"
	with open(os.path.join(args.output_dir, filename), 'w') as dest:
		if s_format == "human":
			for k, value in stats.items():
				print(f"{k}\t{value}", file=dest)
		if s_format == "json":
			print(json.dumps(stats), file=dest)
		if s_format == "table":
			header = "\t".join(stats.keys())
			data = "\t".join(list(map(str, stats.values())))
			print(f"#{header}", file=dest)
			print(f"{data}", file=dest)
	return filename


def assemble(args, reads, estimated_genome_size) -> str:
	'''Returns path to assembled genome'''
	logger = logging.getLogger("main")
	logger.info("Assembling started")
	aslydir = os.path.join(args.output_dir, "assembly")

	if args.assembler == "flye":
		assembly = flye_assemble(args, reads, estimated_genome_size, aslydir)
	elif args.assembler == "spades":
		assembly = spades_assemble(args, reads, aslydir)
	elif args.assembler == "unicycler":
		assembly = unicycler_assemble(args, reads, aslydir)
	elif args.assembler == "megahit":
		assembly = megahit_assemble(args, reads, aslydir)
	else:
		logger.critical("Not yet implemented")
		raise Exception("Unknown option for assembly")

	stats = assembly_stats(assembly)
	logger.info("Assembly length: %s", stats['Total length'])
	logger.info("Contig count: %s", stats['Sequence count'])
	logger.info("N50: %s", stats['N50'])
	logger.info("L50: %s", stats['L50'])
	logger.info("N90: %s", stats['N90'])
	logger.info("L90: %s", stats['L90'])
	write_assembly_stats(args, stats, prefix=args.assembler, s_format="table")

	return assembly
