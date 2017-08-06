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
    sec_num=0
    text_section="\n\n## Each section must have a title, that will form the table of content."
    text_section+="\n## This table of content is automatically generated. A tag is also requested, to allow easier maintenance of external links."
    text_section+="\n## Note the small (one space) indentation for the title and body keys.\n"
    title_tag='\n title: "To be filled"\n tag: \n body: |\n'
    for i_old,line in enumerate(file_old):
      if "intro : |" in line:
        file_new[0]="## This YAML file contains the introduction as well as the body of the html lesson.\n"
        file_new[1]="## In order to modify the other parts, modify the file lessons.html .\n"
        file_new[2]="\n## This is the introduction ...\n"
        file_new.append(line)
        i_new+=1
      elif "This is the body" in line:
        line="\n## Now comes the different sections, numbered."+text_section
        file_new.append(line)
        i_new+=1
      elif "body : |" in line:
        line='sec0:'+title_tag
        file_new.append(line)
        sec_num=1
      elif "name=" in line:
        line=text_section+'sec%s:'%(sec_num)+title_tag+line
        file_new.append(line)
        i_new+=1
        sec_num+=1
      else:
        file_new.append(line)
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
