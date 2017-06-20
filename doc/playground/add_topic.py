#!/usr/bin/env python
# Python and numpy are required to use this script.

from sys import argv, exit
import numpy as np
import os

filnam="abinit_vars.yml"
Input= file(filnam,'r')
outfile = open("abinit_vars_out.yml", "wb+")
#outfile1 = open("varname_in.yml", "wb+")
for oline in Input:
  if oline[0:11] not in "    varname":
    outfile.write(oline)
  if oline[0:11] in "    varname":
    outfile.write("    topic_class:"+"\n")
    outfile.write("    topic_name:"+"\n")
    outfile.write(oline)

outfile.close()




