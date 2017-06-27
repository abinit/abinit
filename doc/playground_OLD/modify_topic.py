#!/usr/bin/env python
# Python and numpy are required to use this script.

from sys import argv, exit
import numpy as np
import os

filnam="varname_in.yml"
filnam1="varname_out.yml"
Input= file(filnam,'r')
outfile = open(filnam1, "wb+")

for oline in Input:
  if oline[0:11] in "    abivarname":
    oline2=Input.next()
    if oline2[0:15] in "    topic_class":
      cl=oline2[17:27]
      cl1=cl.strip()
      oline3=Input.next()
      nam=oline3[16:33]
      if nam=="":
        outfile.write(oline)
        outfile.write("    topics:"+"\n")
      else:
        outfile.write(oline)
        nam1=nam.strip()
        outfile.write("    topics: "+nam1+"_"+cl1+"\n")
    if oline2[0:10] in "    topics":
      outfile.write(oline)
      outfile.write(oline2)


outfile.close()




