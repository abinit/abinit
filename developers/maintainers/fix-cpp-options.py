#!/usr/bin/env python

import os
import re
import sys

renamed_cppopts = {
  "HAVE_ALGO_LEVMAR":"HAVE_LEVMAR",
  "HAVE_DFT_ATOMPAW":"HAVE_ATOMPAW",
  "HAVE_DFT_BIGDFT":"HAVE_BIGDFT",
  "HAVE_DFT_LIBXC":"HAVE_LIBXC",
  "HAVE_DFT_WANNIER90":"HAVE_WANNIER90",
  "HAVE_DFT_WANNIER90_V1":"HAVE_WANNIER90_V1",
  "HAVE_MATH_MLIB":"HAVE_LINALG_MLIB",
  "HAVE_TIMER_PAPI":"HAVE_PAPI",
  "HAVE_TRIO_ETSF_IO":"HAVE_ETSF_IO",
  "HAVE_TRIO_NETCDF":"HAVE_NETCDF",
  "HAVE_TRIO_NETCDF_MPI":"HAVE_NETCDF_MPI",
  "HAVE_TRIO_PSML":"HAVE_PSML",
  "HAVE_TRIO_YAML":"HAVE_YAML"}
renamed_keys = renamed_cppopts.keys()
renamed_keys.sort()

m4 = re.compile("\.m4$")
chdr = re.compile("\.[ch]$")
fortran = re.compile("\.([Ff]|[Ff]90|finc)$")
pynit = re.compile("__init__.py$")
tinp = re.compile("\.in$")
cppline = re.compile("^#")
cppkeys = ("define .*","include.*","ifdef","ifndef","elif","^if ","else","endif","defined","undef","!","&&","\|\|","\(","\)")
m4line = re.compile("AC_DEFINE")
tline = re.compile("need_cpp_vars")

def fix_cppopts(top):
  cppopts = dict()
  for root,dirs,files in os.walk(top):
    for src in files:
      is_fort = False
      if ( fortran.search(src) or chdr.search(src) ):
        is_fort = True
      is_m4 = False
      if ( m4.search(src) ):
        is_m4 = True
      is_pynit = False
      if ( pynit.search(src) ):
        is_pynit = True
      is_tinp = False
      if ( tinp.search(src) ):
        is_tinp = True
      if ( is_fort or is_m4 or is_pynit or is_tinp ):
        code = file(os.path.join(root,src),"r").readlines()
        cmod = False

        for cidx in range(len(code)):
          line = code[cidx]
          lmod = False

          if ( is_fort and cppline.match(line) ):
            line = re.sub("^#","",line).strip()
            for kw in cppkeys:
              line = re.sub(kw,"",line)
            lmod = True
          elif ( is_m4 and m4line.search(line) ):
            line = re.sub(".*AC_DEFINE","", line)
            line = re.sub("AC_DEFINE","", line)
            line = re.sub(",.*","",line)
            line = re.sub("[\(\[\]\)]","",line)
            lmod = True
          elif ( is_pynit ):
            line = re.sub("[ \",]","",line)
            lmod = True
          elif ( is_tinp and tline.search(line) ):
            line = re.sub(".*need_cpp_vars.*=[ ]*","",line)
            lmod = True

          if ( lmod ):
            line = line.split()
            for item in line:
              if ( item in renamed_keys ):
                cmod = True
                code[cidx] = re.sub(item,renamed_cppopts[item],code[cidx])

        if ( cmod ):
          file(os.path.join(root,src),"w").write("".join(code))

  return 0


if __name__ == "__main__":

  if len(sys.argv) == 1: 
    top = "src"
  else:
    top = sys.argv[1] 

  exit_status = fix_cppopts(top)
  sys.exit(exit_status)
