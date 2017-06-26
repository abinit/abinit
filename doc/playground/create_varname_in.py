#!/usr/bin/env python
# Python and numpy are required to use this script.

from sys import argv, exit
import numpy as np
import os

filnam="abinit_vars.yml"
filnam1="varname_in.yml"
Input= file(filnam,'r')
outfile = open(filnam1, "wb+")

for oline in Input:
  if oline[0:14] in "    abivarname":
    outfile.write(oline)
    outfile.write("    topic_class:"+"\n")
    outfile.write("    topic_name:"+"\n")


outfile.close()




