#!/usr/bin/env python
# Python and numpy are required to use this script.

from sys import argv, exit
import numpy as np
import os

filnam="abinit_vars_play.yml"
filnam1="varname_in.yml"
Input= file(filnam,'r')
Input1= file(filnam1,'r')
outfile = open("abinit_vars_out_play.yml", "wb+")

for oline in Input:
  if oline[0:11] not in "    abivarname":
    outfile.write(oline)
  if oline[0:11] in "    abivarname":
    for oline1 in Input1:
      if oline1 == oline:
        oline2=Input1.next()
        outfile.write(oline2)
        outfile.write(oline)
        break
outfile.close()




