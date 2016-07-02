#!/usr/bin/python

import csv
import os
import re

ign_keywords = [
  "elemental",
  "integer",
  "logical",
  "pure",
  "recursive"]

with open("../abinit-orphans-7.9.1-raw.csv", "rb") as csvfile:
  for row in csv.reader(csvfile):
    my_dir, my_src, my_sub, my_typ, my_act, my_sta = row
    if ( my_typ == "" ):
      sub_found = False
      with open(os.path.join("src", my_dir, my_src), "r") as fh:
        lines = fh.readlines()
      for line in lines:
        line = re.sub("!.*", "", line).lower()
        if ( re.search(my_sub.lower(), line) ):
          line = line.strip().split()
          if ( line[0] in ign_keywords ):
            line = line[1:]
          sub_type = line[0]
          try:
            sub_name = re.sub("\\(.*", "", line[1])
          except:
            sub_name= ""
          if ( sub_name == my_sub.lower() ):
            if ( sub_type == "subroutine" ):
              if ( my_sub == re.sub(".F90", "", my_src) ):
                my_act = "quarantine"
              else:
                my_act = "prune"
              my_sta = "pending"
              print("%s,%s,%s,orphan,%s,%s" % \
                (my_dir, my_src, my_sub, my_act, my_sta))
              sub_found = True
              break
            elif ( sub_type == "function" ):
              print("%s,%s,%s,function,skip,done" % \
                (my_dir, my_src, my_sub))
              sub_found = True
              break
      if ( not sub_found ):
        if ( my_sub != re.sub(".F90", "", my_src) ):
          my_act = "prune"
          my_sta = "done"
          print("%s,%s,%s,orphan,%s,%s" % \
            (my_dir, my_src, my_sub, my_act, my_sta))
        else:
          print(",".join(row))
    else:
      print(",".join(row))
