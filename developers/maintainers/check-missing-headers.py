#!/usr/bin/env python
from __future__ import print_function

import os
import re
import sys

re_f90 = re.compile("\.[Ff]90$")
re_ign = re.compile(r"_(ext|noabirule)")
re_int = re.compile("interfaces_")

def main(source_tree):

  tests = dict()
  tests["mod"] = re.compile(r"(module|program)",re.IGNORECASE)
  tests["sub"] = re.compile(r"end (subroutine|function)",re.IGNORECASE)
  tests["hdr"] = re.compile(r"!!\*\*\*\*")
  tests["fun"] = re.compile("!! FUNCTION")

  ign = list()
  sbm = list()

  cnt = dict()
  acc = dict()
  glo = dict()
  dcc = dict()

  glo["mod"] = 0
  glo["sub"] = 0
  glo["hdr"] = 0
  glo["fun"] = 0
  glo["dsc"] = 0
  glo["mis"] = 0

  print( "Status of RoboDOC headers in %s" % (source_tree))
  print( "")
  print( "")
  print( "")

  for root,dirs,files in os.walk(source_tree):
    acc["mod"] = 0
    acc["sub"] = 0
    acc["hdr"] = 0
    acc["fun"] = 0
    acc["dsc"] = 0
    acc["mis"] = 0

    dirs.sort()
    files.sort()
    if ( "01_interfaces_ext" in dirs ):
      ign.append(os.path.join(root,"01_interfaces_ext"))
      dirs.remove("01_interfaces_ext")
    for d in dirs:
      if ( re_ign.search(d) ):
        ign.append(os.path.join(root,d))
        dirs.remove(d)

    sources = list()
    if ( re_ign.search(root) ):
      ign.append(root)
    else:
      for src in files:
        if ( (re_f90.search(src)) and not (re_int.match(src)) ):
          sources.append(os.path.join(root,src))

    if ( len(sources) > 0 ):
      print ("File                                             Mods Subs Hdrs Funs Dscs Miss")
      print ("------------------------------------------------ ---- ---- ---- ---- ---- ----")

    for src in sources:
      cnt["mod"] = 0
      cnt["sub"] = 0
      cnt["hdr"] = 0
      cnt["fun"] = 0
      cnt["dsc"] = 0
      cnt["mis"] = 0

      chk_fun = False
      with open(src, "r") as fh:
        inp = fh.readlines()

      for line in inp:
        tmp = line.split()
        if ( chk_fun ):
          if ( len(tmp) > 5 ):
            cnt["dsc"] += 1
          chk_fun = False

        for blk in tests:
          if ( tests[blk].match(line.strip()) ):
            cnt[blk] += 1
            if ( blk == "fun" ):
              chk_fun = True

      cnt["mis"] = (cnt["mod"]+cnt["sub"]-cnt["fun"]) + (cnt["fun"]-cnt["dsc"])
      if ( cnt["mis"] < 0 ):
        cnt["mis"] = 0

      print ("%-48s %4d %4d %4d %4d %4d %4d" % \
        (src,cnt["mod"],cnt["sub"],cnt["hdr"],cnt["fun"],cnt["dsc"],cnt["mis"]))

      for blk in cnt:
        acc[blk] += cnt[blk]

      if ( (cnt["sub"] > 1) and (cnt["mod"] == 0) ):
        sbm.append(src)

    if ( len(sources) > 0 ):
      my_sum = acc["mod"] + acc["sub"]
      if ( my_sum > 0 ):
        dcc[root] = int(float(my_sum-acc["mis"])*100.0/float(my_sum))

      print ("------------------------------------------------ ---- ---- ---- ---- ---- ----")
      print ("Total in %-39s %4d %4d %4d %4d %4d %4d" % \
        (root,acc["mod"],acc["sub"],acc["hdr"],acc["fun"],acc["dsc"],acc["mis"]))
      print ("")
      print ("")
      print ("")

    for blk in cnt:
      glo[blk] += acc[blk]

  print ("================================================ ==== ==== ==== ==== ==== ====")
  print ("Summary                                          Mods Subs Hdrs Funs Dscs Miss")
  print ("------------------------------------------------ ---- ---- ---- ---- ---- ----")
  print ("%-48s %4d %4d %4d %4d %4d %4d" % \
    ("*",glo["mod"],glo["sub"],glo["hdr"],glo["fun"],glo["dsc"],glo["mis"]))
  print ("================================================ ==== ==== ==== ==== ==== ====")
  print ("")
  print ("Mods: number of Fortran modules or programs")
  print ("Subs: number of subroutines or functions")
  print ("Hdrs: number of RoboDOC headers")
  print ("Funs: number of RoboDOC 'FUNCTION' keywords")
  print ("Dscs: number of actually described functions")
  print ("Miss: number of missing descriptions")
  print ("")
  print ("")
  print ("")

  print ("Directory                  % done")
  print ("------------------------ --------")
  dcc_keys = sorted(dcc.keys())
  for key in dcc_keys:
    print ("%-24s %8d" % (key,dcc[key]))
  my_sum = glo["mod"] + glo["sub"]
  if ( my_sum > 0 ):
    pct = int(float(my_sum-glo["mis"])*100.0/float(my_sum))
    print ("------------------------ --------")
    print ("%-24s %8d" % ("*",pct))
  print ("")
  print ("")
  print ("")

  print ("The following directories have been ignored:")
  print ("")
  for d in ign:
    print ("  * %s/" % (d))
  print ("")
  print ("")
  print ("")

  print ("The following files are not modules, but contain a collection of routines:" )
  print ("")
  for src in sbm:
    print ("  * %s" % (src))
  print ("")
  print ("")
  print ("")

if __name__ == "__main__":

  if len(sys.argv) == 1: 
    source_tree = "src"
  else:
    source_tree = sys.argv[1] 

  main(source_tree)
