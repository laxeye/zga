#!/usr/bin/env python3
import argparse
import os.path
import logging
import shutil
import subprocess
import re

def parse_args():
	parser = argparse.ArgumentParser(description = "ZGA genome annotator")
	
	#General options
	parser.add_argument("-s", "--step", help="Starting step of pipeline", 
		default="qc", choices=["qc", "processing", "assembling", "annotating", "checking"])
	parser.add_argument("-o", "--output-dir", default="zga-output", help="Output directory")
	parser.add_argument("-t", "--threads", type=int, default=1, 
		help="Number of CPU threads to use (where possible)")
	parser.add_argument("-m", "--memory-limit", type=int, default=16, 
		help="Memory limit in GB")
	parser.add_argument("--tmp-dir", default="zga-temp", help="Temporary directory")
	parser.add_argument("--genus", default="Unknown", help="")
	parser.add_argument("--species", default="sp.", help="")
	parser.add_argument("--strain", default="", help="")
	parser.add_argument("--force", action="store_true", help="")
	# parser.add_argument("--", help="")


	#Input
	parser.add_argument("-1", "--pe-1", help="")
	parser.add_argument("-2", "--pe-2", help="")
	parser.add_argument("-M", "--pe-merged", help="")
	parser.add_argument("-S", "--single-end", help="")
	parser.add_argument("--mp-1", help="Mate pair forward reads")
	parser.add_argument("--mp-2", help="Mate pair forward reads")
	# parser.add_argument("--mp-i", help="Mate pair interleaved reads")
	parser.add_argument("--pacbio", help="")
	parser.add_argument("--nanopore", help="")
	
	# Read processing
	parser.add_argument("-q", "--quality-cutoff", type=int, default=25, 
		help="Base quality cutoff")
	parser.add_argument("--adapters", default="util/TruSeq.adapters.fna", 
		help="Adapter sequences for trimming")
	#parser.add_argument("-", "--", help="")
	
	#Assembler
	parser.add_argument("-a", "--assembler", default="spades", choices=["spades", "unicycler"])
	parser.add_argument("--no-correct", action="store_true", help="Disable read correction")

	# Should be implemented
	# parser.add_argument("--interleaved", help="")
	# parser.add_argument("--phred-offset", help="")

	# Annotation
	parser.add_argument("--locus-tag", default="", help="Locus tag")
	parser.add_argument("--gcode", default=11, type=int, help="Genetic code")
	parser.add_argument( "--center_name", default="", help="")

	#parser.add_argument("-", "--", default="", help="")

	#group = parser.add_mutually_exclusive_group(required=True)
	#group.add_argument("-l", "--low_triangular_matrix", help="Low triangular matrix of distance values (ANI or AAI)")
	#group.add_argument("-t", "--table", help="Tab separated table of distance values")
	
	return parser.parse_args()
'''
def missing_file(file):
	logging.error("Missing file: %s!" % file)
'''

def check_reads(args):
	#logger = logging.getLogger("check_reads")
	reads = {}
	read_list = [args.pe_1, args.pe_2, args.single_end, args.pe_merged, 
	args.mp_1, args.mp_2, args.pacbio, args.nanopore]

	logger.info("Checking input files.")

	for f in read_list:
		if f and not os.path.isfile(f):
			logger.error("File %s doesn't exist" % f)
			raise FileNotFoundError("File %s doesn't exist" % f)
			exit(1)

	if args.pe_1:
		if args.pe_2:
			reads["pe_1"] = os.path.abspath(args.pe_1)
			reads["pe_2"] = os.path.abspath(args.pe_2)
		else:
			reads["single"] = os.path.abspath(args.pe_1)
			logger.warning("Only forward reads available")

	if args.pe_merged:
		reads["merged"] = os.path.abspath(args.pe_merged)

	return reads

def create_subdir(subdir):
	path = os.path.join(output_dir, subdir)
	try:
		os.mkdir(path)
	except Exception as e:
		raise e("Impossible to create directory \"%s\"" % path + 
			"Check provided path and permisions")	
	return path

def read_QC(args, reads):
	logger.info("Read quality control started")
	qcoutdir = create_subdir("QC")
	cmd = ["fastqc", "-q", "-t", str(args.threads), "-o", qcoutdir] + list(reads.values())
	logger.debug(cmd)
	rc = subprocess.run(cmd, stdout = subprocess.PIPE, stderr = subprocess.DEVNULL).returncode
	logger.info("Read quality control finished")
	return rc

def pe_read_processing(args, reads):
	logger.info("Read processing started")
	readdir = create_subdir("reads")
	
	#Trim and filter
	MINLEN = 55
	WINDOW = 3
	
	# Paired
	if 'pe_1' in reads.keys() and 'pe_2' in reads.keys():
		logger.debug("Trimming and filtering paired end reads")
		out_pe1 = os.path.join(readdir, "pe_1.fq")
		out_pe2 = os.path.join(readdir, "pe_2.fq")
		cmd = ["fastq-mcf", "-H", "-X", "-q", str(args.quality_cutoff), "-l", str(MINLEN), 
			'-w', str(WINDOW), "-o", out_pe1, "-o", out_pe2, "n/a", reads["pe_1"], reads["pe_2"]]
		logger.debug(' '.join(cmd))
		rc = subprocess.run(cmd, stdout = subprocess.PIPE, stderr = subprocess.DEVNULL).returncode
		if rc != 0:
			logger.error("An error during processing paired-end reads")
		else:
			reads['pe_1'] = out_pe1
			reads['pe_2'] = out_pe2

	# Single
	if 'single' in reads.keys():
		logger.debug("Trimming and filtering single end reads")
		out_single = os.path.join(readdir, "single.fq")
		cmd = ["fastq-mcf", "-H", "-X", "-q", str(args.quality_cutoff), "-l", str(MINLEN), 
			"-o", out_single, "n/a", reads["single"]]
		logger.debug(' '.join(cmd))
		rc = subprocess.run(cmd, stdout = subprocess.PIPE, stderr = subprocess.DEVNULL).returncode
		if rc != 0:
			logger.error("An error during processing single-end reads")
		else:
			reads['single'] = out_single

	#Merging with seqprep
	if 'merged' not in reads.keys() and 'pe_1' in reads.keys() and 'pe_2' in reads.keys():
		nmr1 = os.path.join(readdir, "nm.pe_1.fq.gz")
		nmr2 = os.path.join(readdir, "nm.pe_2.fq.gz")
		merged = os.path.join(readdir, "merged.fq.gz")
		cmd = ['seqprep', '-f', reads['pe_1'], '-r', reads['pe_2'], '-1', nmr1, '-2', nmr2, '-s', merged]
		logger.debug(' '.join(cmd))
		rc = subprocess.run(cmd, stdout = subprocess.PIPE, stderr = subprocess.DEVNULL).returncode
		if rc != 0:
			logger.error("An error during merging paired-end reads")
		else:
			reads['merged'] = merged
			reads['pe_1'] = nmr1
			reads['pe_2'] = nmr2

	logger.info("Read processing finished")
	return rc, reads

def assemble(args, reads):
	logger.info("Assembling started")
	outdir = os.path.join(output_dir, "assembly")
	if args.assembler == 'spades':
		#check for version of spades to use "--merged", "--isolate" etc.
		
		spades_help = str(subprocess.run(['spades.py'], stdout = subprocess.PIPE, stderr = subprocess.PIPE,
			universal_newlines=True).stdout)
		version = re.search(r'(?<=v)\d\S+', spades_help).group(0)
		logger.debug("Spades version %s detected" % version)
		v_major = int(version.split(".")[0])
		v_minor = int(version.split(".")[1])
		
		cmd = ['spades.py', '-o', outdir, '--careful', '-t', str(args.threads), 
			'-m', str(args.memory_limit), '--cov-cutoff', 'auto']
		'''
		I'm not sure to use this mode...
		if v_major > 3 or (v_major >= 3 and v_minor >= 14):
			cmd += ['--isolate']
		else:
			cmd += ['--careful']
		'''
		if 'pe_1' in reads.keys() and 'pe_2' in reads.keys():
			cmd += ['-1', reads['pe_1'], '-2', reads['pe_2']]
		if 'merged' in reads.keys():
			cmd += ['--merged', reads['merged']]
		if 'single' in reads.keys():
			cmd += ['-s', reads['single']]
		if 'nanopore' in reads.keys():
			cmd += ['--nanopore', reads['nanopore']]
		if 'pacbio' in reads.keys():
			cmd += ['--pacbio', reads['pacbio']]
		if args.no_correct == True:
			cmd += ['--only-assembler']
		logger.debug(" ".join(cmd))
		subprocess.run(cmd, stderr = subprocess.PIPE, stdout = subprocess.PIPE)
	else:
		logger.critical("Not yet implemented")
	# TODO unicycler

# checkm taxonomy_wf

# annotation DFAST

def main():
	args = parse_args()
	global output_dir

	output_dir = os.path.abspath(args.output_dir)
	if os.path.isdir(output_dir):
		if args.force:
			shutil.rmtree(output_dir)
		else:
			raise FileExistsError("\nOutput directory \"%s\" already exists.\n" % output_dir + 
				"Use --force to overwrite or provide another path")
	try:
		os.mkdir(output_dir)
	except Exception as e:
		raise e("Imposible to create directory \"%s\"" % output_dir + 
			"Check provided path and permisions")

	global logger
	logger = logging.getLogger("main")
	logger.setLevel(logging.DEBUG)
	formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
	ch = logging.StreamHandler()
	ch.setLevel(logging.INFO)
	ch.setFormatter(formatter)
	logger.addHandler(ch)
	fh = logging.FileHandler(os.path.join(output_dir,"zga.log"))
	fh.setLevel(logging.DEBUG)
	fh.setFormatter(formatter)
	logger.addHandler(fh)

	logger.info("Start")

	steps = {"qc":1, "processing":2, "assembling":3, "annotation":4, "check":5}
	start_step_int=steps[args.step]
	if start_step_int <= 3:
		reads = check_reads(args)
		logger.debug(reads)
	if start_step_int == 1:
		status = read_QC(args, reads)
		logger.debug(status)
	if start_step_int <= 2:
		status, reads = pe_read_processing(args, reads)
		logger.debug(status)
		logger.info(reads)
	if start_step_int <= 3:
		assemble(args, reads)
	if start_step_int <= 4:
		logger.info("Annotating genome")
	if start_step_int <= 5:
		logger.info("Checking genome quality")

if __name__ == '__main__':
	main()