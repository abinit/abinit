#! /usr/bin/env python 
"""
   This script is intended as a template for quick modification of large files, e.g. 
   inserting one more field in abinit_vars.yml
   Modify at your will !

   Feed it with a list of files (at least one file is needed)
   python filemodifier file1 file2 file3 ...
"""

from __future__ import print_function

import sys
import os
import yaml
import re
import string
import argparse

################################################################################

# Generic checks

################################################################################

filelist = sys.argv[1:]
if filelist != [] :
  if filelist[0] == "-h" or filelist[0] == "--help" :
    print(__doc__)
    sys.exit()
else:
  print(__doc__)
  sys.exit()

################################################################################

for path in filelist:
  print(" Working on file '%s'."%path)
  with open(path,'r') as f_old:
    file_old=f_old.readlines()
    f_new=open(path+"_new","w")

################################################################################

    #####################################################
    #The file transformations can be made in this section

    for line in file_old:
      f_new.write(line)
      if "excludes" in line:
        f_new.write("    executables: abinit\n")

    #####################################################

################################################################################

  #Finishes by moving the new file at the new place
  print(" Changes done")
  os.system("mv %s %s_old"%(path,path))
  os.system("mv %s_new %s"%(path,path))
  os.system("rm %s_old"%(path))
  print(" New file %s written"%path)
