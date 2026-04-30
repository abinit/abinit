#!/usr/bin/env python
"""
Check LibPAW bindings.

This script automates the creation of the libpaw-bindings tarball by
running the configuration and make process. It then extracts and builds
the tarball in a temporary directory to verify its integrity, and cleans
up the temporary files afterward.
"""

import glob
import os
import sys
import tempfile
from subprocess import Popen

try:
    sys.setdefaultencoding("utf8")
except AttributeError:
    pass
from shutil import rmtree


def main(home_dir: str = "") -> int:
  """
  Main logic for testing the libpaw tarball creation and build process.

  Args:
      home_dir: The base directory for libpaw bindings (default: "").

  Returns:
      0 if the tarball was created and compiled successfully, 1 otherwise.
  """
  # create tarball
  sys.stdout.write("Creating tarball...\n")

  if home_dir!="":
    bdir = os.path.join(home_dir,"bindings")
    cmd = "cd "+bdir+" && ./configure && make libpaw-bindings && cp libpaw/*tar.gz /tmp"
  else:
    cmd = "cd bindings && ../../bindings/configure && make libpaw-bindings && cp libpaw/*tar.gz /tmp"
  ou = tempfile.TemporaryFile()
  er = tempfile.TemporaryFile()

  process = Popen(cmd, shell=True, stdout=ou, stderr=er)
  process.wait()
  process.communicate()
  rc=process.returncode
  if rc != 0:
     ou.seek(0)
     er.seek(0)
     sys.stdout.write("%s\n" % ou.read())
     sys.stderr.write("%s\n" % er.read())
     ou.close()
     er.close()
     retval = 1
     return retval

  ou.close()
  er.close()
  sys.stdout.write(" done...\n")

  # test tarball
  sys.stdout.write("\nTesting tarball...\n")

  cmd = "cd /tmp;libpaw=`ls libpaw*.tar.gz`; tar xzf $libpaw; cd libpaw; make"
  ou = tempfile.TemporaryFile()
  er = tempfile.TemporaryFile()

  process = Popen(cmd, shell=True, stdout=ou, stderr=er)
  process.wait()
  process.communicate()
  rc=process.returncode
  if rc != 0:
     ou.seek(0)
     er.seek(0)
     sys.stdout.write("%s\n" % ou.read())
     sys.stderr.write("%s\n" % er.read())
     ou.close()
     er.close()
     retval = 1
     return retval

  ou.seek(0)
  er.seek(0)
  sys.stdout.write("%s\n" % ou.read())
# sys.stderr.write("%s\n" % er.read())
  ou.close()
  er.close()

  # cleaning
  retval = 0
  try:
    rmtree("/tmp/libpaw")
    f = glob.glob("/tmp/libpaw*.tar.gz")
    os.remove(f[0])
  except:
    sys.stderr.write("cleaning error")
    retval = 1

  return retval

if __name__ == "__main__":
  if len(sys.argv) == 1:
    home_dir = "."
  else:
    home_dir = sys.argv[1]

  sys.exit(main(home_dir))
