#!/usr/bin/python
from __future__ import print_function

import csv
import os
import sys

# Get directory to check as argument
if ( len(sys.argv) > 1 ):
  rf_dir = sys.argv[1]
  if ( not os.path.isdir(os.path.join("src", rf_dir)) ):
    sys.stderr.write("Error: %s is not a directory\n" % rf_dir)
    sys.exit(1)
else:
  sys.exit(0)

# Init
my_cmds = ""
my_orphans = []

# Find source files to quarantine
with open("../abinit-orphans-7.9.1.csv", "rb") as csvfile:
  for row in csv.reader(csvfile):
    my_dir, my_src, my_sub, my_typ, my_act, my_sta = row
    if ( my_dir == rf_dir ):
      if ( (my_typ == "orphan") and \
           (my_act == "quarantine") and \
           (my_sta == "pending") ):
        my_orphans.append(my_src)
        my_cmds += "Remember to remove %s from the src/%s/abinit.src file\n" % \
          (my_src, my_dir)

# Quarantine files
my_orphans = list(set(my_orphans))
for my_src in my_orphans:
  os.system("bzr mv src/%s/%s src/quarantine/%s_%s" % \
    (rf_dir, my_src, rf_dir, my_src))

# Display instructions
my_cmds += "Remember to mark quarantined routines as done\n"
my_cmds += "bzr commit -m \"Quarantined orphan routines from %s\"\n" % rf_dir
print(my_cmds)
