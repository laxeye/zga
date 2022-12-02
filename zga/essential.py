import logging
import os.path
import subprocess
import sys


def create_subdir(parent, child) -> str:
	'''Tries to create a subdirectory, returns path if succes.'''
	logger = logging.getLogger("main")
	path = os.path.join(parent, child)
	try:
		os.mkdir(path)
	except Exception as ex:
		logger.critical("Impossible to create directory \"%s\"", path)
		raise ex
	return path


def run_external(args, cmd, keep_stdout=False, keep_stderr=True):
	'''Run external command using subprocess.

	Returns subprocess.CompletedProcess
	'''
	logger = logging.getLogger("main")
	logger.debug("Running: %s", " ".join(cmd))
	stderr_dest = subprocess.PIPE if keep_stderr else subprocess.DEVNULL
	stdout_dest = subprocess.PIPE if keep_stdout else subprocess.DEVNULL

	try:
		compl_proc = subprocess.run(cmd,
			check=True, stderr=stderr_dest, stdout=stdout_dest, encoding="utf-8")
		if args.transparent:
			print(compl_proc.stderr, file=sys.stderr)
		return compl_proc
	except subprocess.CalledProcessError as err:
		logger.error("Error during execution of: %s", ' '.join(err.cmd))
		logger.info("Please see the logfile for additional information: %s",
			os.path.join(args.output_dir, "zga.log"))
		logger.debug("External tool stderr:\n%s", err.stderr)
		return None
