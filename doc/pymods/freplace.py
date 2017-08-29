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

# Generic checks
filelist = sys.argv[1:]
if filelist != [] :
  if filelist[0] == "-h" or filelist[0] == "--help" :
    print(__doc__)
    sys.exit()
else:
  print(__doc__)
  sys.exit()

# Loop on files
for path in filelist:
  print(" Working on file '%s'."%path)

  # Transfer the content of the "old" file to a list of lines, file_old
  with open(path,'r') as f_old:
    file_old=f_old.readlines()

    #####################################################
    #The file transformations can be made in this section, working directly
    #on the list of lines ...

    file_new=[]
    i_new=0
    for i_old,line in enumerate(file_old):
      if 1 :
        line="  "+line
        #line=line.replace("anavarname","abivarname").rstrip()+"@anaddb\n"
      file_new.append("  "+line)
      i_new+=1
    #####################################################

  #Open the new file, and write the content of file_new
  f_new=open(path+"_new","w")
  for line in file_new:
     f_new.write(line)

  #Finishes by possibly moving the new file at the new place
  print(" Changes done")
  if 0:
    os.system("mv %s %s_old"%(path,path))
    os.system("mv %s_new %s"%(path,path))
    os.system("rm %s_old"%(path))
    print(" New file %s written"%path)
  else:
    print(" New file %s_new written"%path)
