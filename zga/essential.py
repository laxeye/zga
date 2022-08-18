import logging
import os.path
import subprocess
import sys


def run_external(args, cmd, keep_stdout=False, keep_stderr=True):
	'''Run external command using subprocess.

	Returns subprocess.CompletedProcess
	'''

	logger = logging.getLogger("main")
	logger.debug("Running: %s", " ".join(cmd))
	stderr_dest = subprocess.PIPE if keep_stderr else subprocess.DEVNULL
	stdout_dest = subprocess.PIPE if keep_stdout else subprocess.DEVNULL

	try:
		r = subprocess.run(cmd,
			check=True, stderr=stderr_dest, stdout=stdout_dest, encoding="utf-8")
		if args.transparent:
			print(r.stderr, file=sys.stderr)
		return r
	except subprocess.CalledProcessError as e:
		logger.error("Error during execution of: %s", ' '.join(e.cmd))
		logger.info("Please see the logfile for additional information: %s",
			os.path.join(args.output_dir, "zga.log"))
		logger.debug("External tool stderr:\n%s", e.stderr)
		return None
