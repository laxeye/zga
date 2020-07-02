#!/usr/bin/env python3
import argparse
import os.path
import logging
import shutil
import subprocess
import re
from Bio import SeqIO
import hashlib


def parse_args():
	parser = argparse.ArgumentParser(description="ZGA genome assembly and annotation pipeline")

	# General options
	general_args = parser.add_argument_group(title="General options", description="")
	general_args.add_argument("-s", "--first-step", help="First step of the pipeline", default="qc",
		choices=["qc", "processing", "assembling", "check_genome", "annotation"])
	general_args.add_argument("-l", "--last-step", help="Last step of the pipeline", default="annotation",
		choices=["qc", "processing", "assembling", "check_genome", "annotation"])
	general_args.add_argument("-o", "--output-dir", required=True, help="Output directory")
	general_args.add_argument("--force", action="store_true", help="Overwrite output directory if exists")
	# parser.add_argument("--tmp-dir", default="zga-temp", help="Temporary directory")
	general_args.add_argument("-t", "--threads", type=int, default=1,
		help="Number of CPU threads to use (where possible)")
	general_args.add_argument("-m", "--memory-limit", type=int, help="Memory limit in GB for SPAdes")
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
	# parser.add_argument("--phred-offset", help="")
	input_args.add_argument("--pacbio", help="PacBio reads")
	input_args.add_argument("--nanopore", help="Nanopore reads")

	# Read processing
	reads_args = parser.add_argument_group(title="Setting for processing of reads")
	reads_args.add_argument("-q", "--quality-cutoff", type=int, default=25,
		help="Base quality cutoff for short reads, default: 25")
	reads_args.add_argument("--adapters", 
		help="FASTA file with adapter sequences for trimming from short reads. By default Illumina adapter sequences are used.")
	reads_args.add_argument("--merge-with", default="bbmerge", choices=["bbmerge", "seqprep"],
		help="Tool for merging overlapping paired-end reads: bbmerge (default) or seqprep")
	reads_args.add_argument("--filter-by-tile", action="store_true",
		help="Filter reads based on positional quality over a flowcell.")
	reads_args.add_argument("--genome-size-estimation", action="store_true",
		help="Estimate genome size with mash")
	#Mate pair read processing
	reads_args.add_argument("--use-unknown-mp", action="store_true",
		help="Include reads that are probably mate pairs (default: only known MP used)")

	# Assembly
	asly_args = parser.add_argument_group(title="Assembly settings")
	asly_args.add_argument("-a", "--assembler", default="unicycler", choices=["spades", "unicycler"],
		help="Assembler: unicycler (default; better quality, may use only long reads) or spades (faster, may use mate-pair reads).")
	asly_args.add_argument("--no-correction", action="store_true", help="Disable read correction")
	# Spades options
	asly_args.add_argument("--use-scaffolds", action="store_true",
		help="SPAdes: Use assembled scaffolds.")
	asly_args.add_argument("--spades-k-list",
		help="SPAdes: List of kmers, comma-separated even numbers e.g. '21,33,55,77'")
	# Unicycler options
	asly_args.add_argument("--unicycler-mode", default="normal", choices=['conservative', 'normal', 'bold'],
		help="Unicycler: assember mode: conservative, normal (default) or bold.")
	asly_args.add_argument("--linear-seqs", default=0, help="Expected number of linear sequences")
	asly_args.add_argument("--extract-replicons", action="store_true",
		help="Unicycler: extract replicons (e.g. plasmids) from the assembly to separate files")

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
	anno_args.add_argument("--minimum-length", help="Minimum sequence length in genome assembly.")

	return parser.parse_args()


def check_reads(args):
	reads = {}
	read_list = [args.pe_1, args.pe_2, args.single_end, args.pe_merged,
	args.mp_1, args.mp_2, args.pacbio, args.nanopore]

	logger.info("Checking input files.")

	for f in read_list:
		if f and not os.path.isfile(f):
			logger.error("File %s doesn't exist" % f)
			raise FileNotFoundError("File %s doesn't exist" % f)

	if args.pe_1 and args.pe_2:
		reads['pe_1'] = os.path.abspath(args.pe_1)
		reads['pe_2'] = os.path.abspath(args.pe_2)
	elif args.pe_1 or args.pe_2:
		logger.error("Single end reads provided as paired. Please use \"--single-end\" option")
		exit(1)

	if args.mp_1 and args.mp_2:
		reads['mp_1'] = os.path.abspath(args.mp_1)
		reads['mp_2'] = os.path.abspath(args.mp_2)

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
		raise e(f"Impossible to create directory: {path}")
	return path


def run_external(args, cmd, return_code=True):
	logger.debug("Running: " + " ".join(cmd))
	if not return_code or not args.transparent:
		r = subprocess.run(cmd, stderr=subprocess.PIPE, stdout=subprocess.PIPE, encoding="utf-8")
	else:
		r = subprocess.run(cmd)
	if r.returncode != 0:
		logger.error(f'Non-zero return code of "{" ".join(cmd)}"')

	if return_code:
		return r.returncode
	else:
		return r


def read_QC(args, reads):
	logger.info("Read quality control started")
	qcoutdir = create_subdir(args.output_dir, "QC")
	cmd = ["fastqc", "-q", "-t", str(args.threads), "-o", qcoutdir] + list(reads.values())

	return run_external(args, cmd)


def filter_by_tile(args, reads, readdir):
	filtered_pe_r1 = os.path.join(readdir, "filtered_pe_r1.fq.gz")
	filtered_pe_r2 = os.path.join(readdir, "filtered_pe_r2.fq.gz")

	cmd = ["filterbytile.sh", f"in={reads['pe_1']}", f"in2={reads['pe_2']}",
	f"out={filtered_pe_r1}", f"out2={filtered_pe_r2}"]

	rc = run_external(args, cmd)

	if rc == 0:
		reads['pe_1'] = filtered_pe_r1
		reads['pe_2'] = filtered_pe_r2

	return reads


def merge_seqprep(args, reads, readdir):
	notmerged_r1 = os.path.join(readdir, "nm.pe_1.fq.gz")
	notmerged_r2 = os.path.join(readdir, "nm.pe_2.fq.gz")
	merged = os.path.join(readdir, "merged.fq.gz")
	cmd = ["seqprep", "-f", reads['pe_1'], "-r", reads['pe_2'], "-1", notmerged_r1,
		"-2", notmerged_r2, "-s", merged]
	logger.info("Merging paired-end reads.")

	rc = run_external(args, cmd)

	if rc == 0:
		reads['merged'] = merged
		reads['pe_1'] = notmerged_r1
		reads['pe_2'] = notmerged_r2

	return reads


def merge_bb(args, reads, readdir):
	# Worth to be args?
	bb_trim = True
	bb_trimq = "10"  # should be str

	notmerged_r1 = os.path.join(readdir, "nm.pe_1.fq.gz")
	notmerged_r2 = os.path.join(readdir, "nm.pe_2.fq.gz")
	merged = os.path.join(readdir, "merged.fq.gz")
	cmd = ["bbmerge.sh", f"in1={reads['pe_1']}", f"in2={reads['pe_2']}",
		f"outu1={notmerged_r1}", f"outu2={notmerged_r2}", f"out={merged}"]
	if bb_trim:
		cmd += ["qtrim2=t", f"trimq={bb_trimq}"]
	logger.info("Merging paired-end reads.")

	rc = run_external(args, cmd)

	if rc == 0:
		reads['merged'] = merged
		reads['pe_1'] = notmerged_r1
		reads['pe_2'] = notmerged_r2

	return reads


def trim_and_filter_pe(args, reads, readdir):
	# Trimming and filtering constants
	MINLEN = 33
	WINDOW = 3

	if "pe_1" in reads.keys() and "pe_2" in reads.keys():
		logger.info("Trimming and filtering paired end reads")
		out_pe1 = os.path.join(readdir, "pe_1.fq")
		out_pe2 = os.path.join(readdir, "pe_2.fq")
		cmd = ["fastq-mcf", "-H", "-X", "-q", str(args.quality_cutoff), "-l", str(MINLEN),
			"-w", str(WINDOW), "-o", out_pe1, "-o", out_pe2, args.adapters, reads["pe_1"], reads["pe_2"]]
		rc = run_external(args, cmd)
		if rc == 0:
			reads['pe_1'] = out_pe1
			reads['pe_2'] = out_pe2

	if "single" in reads.keys():
		logger.info("Trimming and filtering single end reads")
		out_single = os.path.join(readdir, "single.fq")
		cmd = ["fastq-mcf", "-H", "-X", "-q", str(args.quality_cutoff), "-l", str(MINLEN),
			"-w", str(WINDOW), "-o", out_single, "n/a", reads["single"]]

		rc = run_external(args, cmd)
		if rc == 0:
			reads['single'] = out_single

	if "merged" in reads.keys():
		logger.info("Trimming and filtering merged paired-end reads")
		out = os.path.join(readdir, "merged.fq")
		cmd = ["fastq-mcf", "-H", "-X", "-q", str(args.quality_cutoff), "-l", str(MINLEN),
			"-w", str(WINDOW), "-o", out, "n/a", reads["merged"]]

		rc = run_external(args, cmd)
		if rc == 0:
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

	if args.filter_by_tile and "pe_1" in reads.keys() and "pe_2" in reads.keys():
		reads = filter_by_tile(args, reads, readdir)

	reads = trim_and_filter_pe(args, reads, readdir)

	# Merging overlapping paired-end reads
	if "merged" not in reads.keys() and "pe_1" in reads.keys() and "pe_2" in reads.keys():
		if args.merge_with == "bbmerge":
			reads = merge_bb(args, reads, readdir)
		else:
			reads = merge_seqprep(args, reads, readdir)

	if "mp_1" in reads.keys():
		reads = mp_read_processing(args, reads, readdir)

	logger.info("Read processing finished")

	return reads


def mash_estimate(args, reads):
	# Minimum copy number of k-mer to include it
	MINIMUM_COPIES = 3

	reads_to_sketch=[]
	short_reads=['pe_1', 'pe_2', 'merged', 'single']
	for key in short_reads:
		if key in reads.keys():
			reads_to_sketch.append(reads[key])
			readdirname=os.path.dirname(reads[key])
	if len(reads_to_sketch) == 0:
		logger.error("Not possible to estimate gemome size: short reads missing")
		return None

	sketchprefix=os.path.join(readdirname, "sketch")	
	cmd = ["mash", "sketch", "-r", "-m", str(MINIMUM_COPIES), "-o", sketchprefix]
	cmd += reads_to_sketch
	r = run_external(args, cmd, False)

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
		logger.info(f"Estimated genome size is {int(best_estimation[0])} at coverage {best_estimation[1]}.")
	else:
		logger.error("Genome size estimation with \"mash\" failed.")

	return r.returncode


def mp_read_processing(args, reads,readdir):
	prefix = os.path.join(readdir, "nxtrim")
	MINLENGTH = 31

	cmd = ["nxtrim", "-1", reads['mp_1'], "-2", reads['mp_2'], "--separate"]
	cmd += ["--justmp", "-O", prefix, "-l", str(MINLENGTH)]
	rc = run_external(args, cmd)
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
		reads['mp_1'] = f"{prefix}_R1.all.fastq.gz"
		reads['mp_2'] = f"{prefix}_R2.all.fastq.gz"
	else:
		reads['mp_1'] = f"{prefix}_R1.mp.fastq.gz"
		reads['mp_2'] = f"{prefix}_R2.mp.fastq.gz"

	return reads


def assemble(args, reads):
	logger.info("Assembling started")
	aslydir = os.path.join(args.output_dir, "assembly")

	if args.assembler == "spades":

		cmd = ["spades.py", "-o", aslydir, "--careful", "-t", str(args.threads),
			"--cov-cutoff", "auto"]

		# Check Spades and get it's version to use verion-specific features
		# as "--merged", "--isolate" etc.
		try:
			spades_version = str(subprocess.run(["spades.py", "-v"], stdout=subprocess.PIPE,
				stderr=subprocess.PIPE, universal_newlines=True).stdout)
		except Exception as e:
			logger.error("Failed to run \"spades\"")
			raise e
		#version = re.search(r'(?<=v)\d\S+', spades_help).group(0)
		version = re.search(r'[\d\.]+', spades_version)[0]
		logger.debug("Spades version %s detected" % version)

		'''
		v_major = int(version.split(".")[0])
		v_minor = int(version.split(".")[1])

		I'm not sure is it worth to use this mode...
		if v_major > 3 or (v_major >= 3 and v_minor >= 14):
			cmd += ['--isolate']
		else:
			cmd += ['--careful']
		'''

		if args.memory_limit:
			cmd += ["-m", str(args.memory_limit)]
		if 'pe_1' in reads.keys() and 'pe_2' in reads.keys():
			cmd += ["-1", reads['pe_1'], "-2", reads['pe_2']]
		if 'merged' in reads.keys():
			cmd += ["--merged", reads['merged']]
		if 'single' in reads.keys():
			cmd += ["-s", reads['single']]
		if 'mp_1' in reads.keys() and 'mp_2' in reads.keys():
			cmd += ["--mp1-1", reads['mp_1'], "--mp1-2", reads['mp_2']]
		if 'nanopore' in reads.keys():
			cmd += ["--nanopore", reads['nanopore']]
		if 'pacbio' in reads.keys():
			cmd += ["--pacbio", reads['pacbio']]
		if args.no_correction:
			cmd += ["--only-assembler"]
		if args.spades_k_list:
			cmd += ["-k", args.spades_k_list]

		rc = run_external(args, cmd)

		if rc != 0:
			logger.error("Genome assembly finished with errors.")
			logger.error("Plese check %s for more information." % os.path.join(aslydir, "spades.log"))
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
			version = subprocess.run(["unicycler", "--version"],
			stderr=subprocess.PIPE, stdout=subprocess.PIPE).stdout
		except Exception as e:
			logger.critical("Failed to run \'Unicycler\'.")
			raise e

		cmd = ["unicycler", "-o", aslydir, "-t", str(args.threads),
			"--mode", args.unicycler_mode]
		if 'pe_1' in reads.keys() and 'pe_2' in reads.keys():
			cmd += ["-1", reads['pe_1'], "-2", reads['pe_2']]
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

		rc = run_external(args, cmd)

		if rc != 0:
			logger.error("Genome assembly finished with errors.")
			logger.error("Plese check %s for more information." % os.path.join(aslydir, "unicycler.log"))
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
	logfile = os.path.join(aslydir, "unicycler.log")
	assemblyfile = os.path.join(aslydir, "assembly.fasta")
	with open(logfile, "r") as log:
	    regexp = re.compile(r'\s?\d+.+\scomplete')
	    repl_lengths = []
	    for l in log.readlines():
	        if regexp.search(l):
	            component, segments, links, length, _, _, status = l.split()
	            if int(segments) == 1:
	                repl_lengths.append(int(length.replace(",","")))
	    logger.debug("Extracting " + str(len(repl_lengths)) + " replicon(s).")
	    with open(assemblyfile, "r") as assembly:
	        replicons = [x for x in SeqIO.parse(assembly, "fasta") if len(x) in repl_lengths]
	        for x in range(len(replicons)):
	            try:
	                F = open(os.path.join(aslydir,f"replicon.{x+1}.fasta"),"w")
	                SeqIO.write(replicons[x], F, "fasta")
	                F.close()
	            except Exception as e:
	                raise e


def locus_tag_gen(genome):
	logger.info("No locus tag provided. Generating it as MD5 hash of genome")
	with open(genome, 'rb') as genomefile:
		digest = hashlib.md5(genomefile.read()).hexdigest()
		locus_tag = "".join([chr(65 + (int(digest[x],16) + int(digest[x+1],16)) % 26) for x in range(0,12,2)])
		logger.info("Locus tag generated: %s" % locus_tag)
		return locus_tag


def annotate(args):
	logger.info("Genome annotation started")

	try:
		version = subprocess.run(["dfast", "--version"], stderr=subprocess.PIPE, stdout=subprocess.PIPE).stdout
	except Exception as e:
		logger.critical("Failed to run DFAST")
		raise e

	annodir = os.path.join(args.output_dir, "annotation")

	cmd = ["dfast", "-g", args.genome, "-o", annodir, "--organism", " ".join([args.genus, args.species]),
		"--cpu", str(args.threads)]
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
	if args.minimum_length:
		cmd += ["--minimum_length", args.minimum_length]

	rc = run_external(args, cmd)

	return os.path.join(annodir, "genome.fna")


def check_phix(args):
	logger.info("Checking assembly for presence of Illumina phiX control")
	phix_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),"data/phiX174.fasta")
	blast_format = "6 sseqid pident slen length"
	cmd = ["blastn", "-query", phix_path, "-subject", args.genome, "-outfmt", blast_format, "-evalue", "1e-6"]

	logger.debug("Running: " + " ".join(cmd))
	try:
		blast_out = str(subprocess.run(cmd, universal_newlines=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE).stdout).rstrip()
	except Exception as e:
		raise e

	phix_contigs = []
	if bool(blast_out):
		for l in blast_out.split('\n'):
			i, p, s, l = l.split("\t")
			if float(p) > 95.0 and int(s)/int(l) > 0.5:
				phix_contigs.append(i)
		phix_contigs = list(set(phix_contigs))

		if len(phix_contigs) > 0:
			logger.info(f"PhiX was found in: {', '.join(phix_contigs)}")
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
		raise e

	checkm_outdir = os.path.join(args.output_dir, "checkm")
	checkm_outfile = os.path.join(args.output_dir, "checkm.result.txt")
	checkm_ext = os.path.splitext(args.genome)[1]

	if args.checkm_mode == "taxonomy_wf":

		cmd = ["checkm", "taxonomy_wf", "-f", checkm_outfile, "-x", checkm_ext, "-t", str(args.threads)]

		try:
			checkm_taxon_list = str(subprocess.run(["checkm", "taxon_list"],
				universal_newlines=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE).stdout)
		except Exception as e:
			raise e

		if args.checkm_taxon and args.checkm_rank:
			found = False
			for l in checkm_taxon_list.split("\n"):
				if args.checkm_taxon in l and args.checkm_rank in l:
					found = True
					break
			if not found:
				logger.error(f'Taxon {args.checkm_taxon} of rank {args.checkm_rank} not available for CheckM')
				args.checkm_taxon = None
				args.checkm_rank = None

		if not args.checkm_taxon or not args.checkm_rank:
			args.checkm_taxon = "Archaea" if args.domain == "archaea" else "Bacteria"
			args.checkm_rank = "domain"

		logger.info(f'{args.checkm_taxon} marker set will be used for CheckM')

		cmd += [args.checkm_rank, args.checkm_taxon, checkm_indir, checkm_outdir]

	else:
		cmd = ["checkm", "lineage_wf", "-f", checkm_outfile, "-t", str(args.threads), "-x", checkm_ext]
		if not args.checkm_full_tree:
			cmd += ["--reduced_tree"]
		cmd += ["--pplacer_threads", str(args.threads), checkm_indir, checkm_outdir]

	rc = run_external(args, cmd)

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
		exit(0)


def main():
	args = parse_args()

	args.output_dir = os.path.abspath(args.output_dir)
	if os.path.isdir(args.output_dir):
		if args.force:
			shutil.rmtree(args.output_dir)
		else:
			raise FileExistsError("\nOutput directory \"%s\" already exists.\n" % args.output_dir +
				"Use --force to overwrite or provide another path")
	try:
		os.mkdir(args.output_dir)
	except Exception as e:
		raise e("Imposible to create directory \"%s\"" % args.output_dir +
			"Check provided path and permisions")

	global logger
	logger = logging.getLogger("main")
	logger.setLevel(logging.DEBUG)
	formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
	ch = logging.StreamHandler()
	ch.setLevel(logging.INFO)
	ch.setFormatter(formatter)
	logger.addHandler(ch)
	fh = logging.FileHandler(os.path.join(args.output_dir,"zga.log"))
	fh.setLevel(logging.DEBUG)
	fh.setFormatter(formatter)
	logger.addHandler(fh)

	logger.info("Pipeline started")

	steps = {"qc":1, "processing":2, "assembling":3, "annotation":5, "check_genome":4}
	args.first_step = steps[args.first_step]
	args.last_step = steps[args.last_step]
	if args.first_step <= 3:
		reads = check_reads(args)
		logger.debug("Reads: " + str(reads))
		if len(list(reads)) == 0:
			logger.error("No reads provided for genome assembly")
			raise Exception("No reads provided for genome assembly")

	if args.first_step > 3 and (not args.genome or not os.path.isfile(args.genome)):
		logger.error("Genome assembly is not provided")
		raise FileNotFoundError()

	# QC
	if args.first_step == 1:
		return_code = read_QC(args, reads)
		check_last_step(args, 1)

	# Processing
	if args.first_step <= 2:
		reads = read_processing(args, reads)
		logger.debug("Processed reads: " + str(reads))
		if args.genome_size_estimation:
			mash_estimate(args, reads)
		check_last_step(args, 2)

	# Assembly
	if args.first_step <= 3:
		args.genome = assemble(args, reads)
		check_last_step(args, 3)

	# Assembly QC
	if args.first_step <= 4:
		logger.info("Checking genome quality")
		if args.check_phix:
			args.genome = check_phix(args)
		checkm_outfile = run_checkm(args)
		if checkm_outfile:
			with open(checkm_outfile) as result:
				completeness, contamination, heterogeneity =  list(map(float, result.readlines()[3].split()[-3::1]))
				if float(completeness) < 80.0 or (float(completeness) - 5.0 * float(contamination) < 50.0):
					logger.info("The genome assembly has low quality!")
				logger.info(f"Genome completeness: {completeness}%")
				logger.info(f"Genome contamination: {contamination}%")
				logger.info(f"Genome heterogeneity: {heterogeneity}%")
		check_last_step(args, 4)

	# Annotation
	if args.first_step <= 5:
		if not args.genome or not os.path.isfile(args.genome):
			logger.error("Genome assembly is not provided")
			raise FileNotFoundError()
		args.genome = annotate(args)
		check_last_step(args, 5)


if __name__ == "__main__":
	main()
