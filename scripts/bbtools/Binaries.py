#!/bin/env python
from __future__ import print_function
import os
import sys
import time
import tempfile
from subprocess import Popen
import configparser

confdir = "config/specs"

binaries = configparser.ConfigParser()
binaries.read(os.path.join(confdir,"binaries.conf"))

general_rc = 0

#for binary in ['newsp']:
for binary in binaries.sections():
	print("**************************************")
	print("Prog ", binary,  ": starting...")
	print("**************************************")
	fcflags = os.environ.get("FCFLAGS", None)
	if fcflags is not None:
		print("WARNING: FCFLAGS is defined in env with value:", fcflags)
		print("Optimization level should be set to -O0 via ac file to speedup compilation!")
	sys.stdout.flush()
	start = time.time()
	cmd = "make -j 10 " + binary
	ou = tempfile.TemporaryFile()
	er = tempfile.TemporaryFile()
	process = Popen(cmd, shell=True, stdout=ou, stderr=er )
	process.wait()
	process.communicate()
	print("\n+++++++  RETURN CODE:", str(process.returncode), "\n")
	print("Completed in:", time.time() - start, "[s]")
	sys.stdout.flush()
	if process.returncode != 0:
		general_rc += 1
		ou.seek(0)
		er.seek(0)
		print(" stdout:\n")
		print(ou.read())
		print(" stderr:\n")
		print(er.read())
		sys.stdout.flush()
	ou.close()
	er.close()

	print("Current Working directory:", os.getcwd())
	filepath = os.path.join(os.getcwd(), "src", "98_main", binary)
	if not os.path.isfile(filepath):
		raise RuntimeError("Cannot find executable:", filepath)

	#cmd = "cd src && make clean > /dev/null 2>&1 && cd .."
	# Run `make clean` in top-level directory to include srd files in shared as well.
	cmd = "make clean > /dev/null 2>&1"
	print(" ----> cleaning process starting...")
	print("command: ", cmd)
	process = Popen(cmd, shell=True )
	process.wait()
	if process.returncode != 0:
		raise RuntimeError("Error while running `%s`:" % cmd)
	print(" ----> command terminated with retcode:", process.returncode)
	print("**************************************\n")
	
	sys.stdout.flush()

sys.exit(general_rc)
