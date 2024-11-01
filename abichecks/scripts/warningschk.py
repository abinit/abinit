#!/usr/bin/env python
from __future__ import division, print_function, absolute_import #unicode_literals, 

import sys
import os
#import re
import os.path
import glob

pack_dir, x = os.path.split(os.path.abspath(__file__))
pack_dir, x = os.path.split(pack_dir)
sys.path.insert(0,pack_dir)
pack_dir, x = os.path.split(pack_dir)
sys.path.insert(0,pack_dir)

from tests.pymods.termcolor import cprint

gnu_warnings = { # ( warning_string, warno, src_excluded )
    #3  : ( 'Unused variable', ['12_hide_mpi','64_psp','68_dmft'] ),
    3  : ( 'Unused variable', [] ),
    4  : ( 'Unused dummy argument',  [] ),
    5  : ( 'Nonstandard type declaration',  ['interfaces','28_numeric_noabirule','01_macroavnew_ext','01_linalg_ext','11_memory_mpi'] ),
    6  : ( 'Same actual argument associated with INTENT', []),  
    7  : ( 'CHARACTER expression will be truncated in assignment',  [] ),
    8  : ( 'Limit of 39 continuations exceeded',  [] ),
    9  : ( 'DOUBLE COMPLEX at (1) does not conform to the Fortran 95 standard',  ['interfaces','01_linalg_ext'] ),
    10 : ( 'at (1) defined but not used', [] ),
    11 : ( 'Character length of actual argument shorter than of dummy argument', [] ),
    #12 : ( 'may be used uninitialized',  [] ), FIXME Disabled cause it segfaults
    13 : ( 'Obsolescent', [] ),
    14 : ( 'Type specified for intrinsic function', [] ),
    15 : ( 'Nonconforming tab character', [] ),
    20 : ( 'Wunused-value', [] ),
}

def usage():
    print("\n Usage: warningschk test_number \n ")


def main(warno, home_dir=""):
  debug = 0

  if not home_dir:
    cwd_dir = os.getcwd()
    if os.path.isabs(sys.argv[0]):
        home_dir = os.path.normpath(os.path.join(os.path.dirname(sys.argv[0]), "../.."))
        inp_dir = os.path.join(home_dir, "abichecks/abirules/Input")
    else:
        inp_dir = os.path.join("..", "Input")
        home_dir = os.path.join(cwd_dir,"../../..")

  else:
    inp_dir = os.path.join(home_dir, "abichecks", "abirules", "Input")
  
  assert os.path.isdir(inp_dir)
  warno = int(warno)
  Warning      = gnu_warnings[warno][0]
  Warning_len  = len(Warning.split(" "))
  src_excluded = gnu_warnings[warno][1]

  # header
  print( "**********************************************************************")
  print( "Warning pattern : '"+Warning+"'")
  print( "**********************************************************************")

  makelog = os.path.join(home_dir, "make.log")
  if not os.path.exists(makelog):
      raise RuntimeError("Cannot find `make.log` file in `%s`.\nUse `make -O multi -j8 > make.log 2>&1`" % home_dir)
  # make.log contains utf-8 characters
  #import io
  #logfile = io.open(makelog, "r", encoding="utf-8")
  logfile = open(makelog)

  words = []
  Buffer = []
  linec = 0
  warning_count = 0
  start = False
  for line in logfile:
      linec = linec + 1
      if linec > 5 : Buffer.pop(0)
      Buffer.append(line)
      if start == False :
          # Examine the make.log file, starting with the section where the directory 10_defs was treated.
          if line.find("Making all in 10_defs") == -1 :
              continue
          else:
              start = True
      if line.find(Warning) != -1 :
          if debug:
              print("[DEBUG] Buffer[0]:", Buffer[0])  # source.F90:line.pos:
              print("[DEBUG] Buffer[2]:", Buffer[2])  # instruction
              print("[DEBUG] Buffer[1]:", Buffer[1])  # position
              print("[DEBUG] Buffer[4]:", Buffer[4])  # Warning: msg
          if True:
              if debug: print("[DEBUG] len of Buffer[0]:", len(Buffer[0].strip()))
              if len(Buffer[0].strip()) != 0:
                  source = Buffer[0].split(":")[0]
                  if source.find('Included at'): source = source.split(" ")[-1]
                  sourceline = Buffer[0].split(":")[1]
                  try:
                      sourceline = sourceline.split(".")[0]
                  except IndexError:
                      pass
                  pattern = os.path.join(home_dir, "src") + "/*/"+source
                  path = glob.glob(pattern)
                  assert len(path) < 2
                  try:
                      source_dir = path[0].split('/')
                      if debug: print ("[DEBUG] source_dir :" + source_dir[-2])
                      if src_excluded.index(source_dir[-2]) :
                          pass
                  except IndexError:
                      pass
                  except ValueError:
                      warning_count += 1
                      try:
                          if warno in [3,4]:
                             warn_msg=Buffer[4].split(" ")[Warning_len+1]
                             print(source + ' = line: ' + sourceline + ', var: ' + warn_msg +' ['+source_dir[-2]+']')
                          elif warno in [6,10]:
                             warn_msg=Buffer[4].split(":")[1].rstrip()
                             warn_code=Buffer[2].rstrip()
                             warn_pos=Buffer[3].rstrip()
                             print("%s = line: %s, " % (source,sourceline),end='')
                             cprint("warn: %s" % (warn_msg),"red")
                             cprint("  ->%s\n  ->%s" % (warn_code,warn_pos),"red")
                          elif warno in [7]:
                             warn_code=Buffer[2].rstrip().lstrip()
                             print("%s = line: %s, " % (source,sourceline),end='')
                             cprint("code: %s" % (warn_code),"red")
                          elif warno in [20]:
                             a = Buffer[4].split(":")[1].split(" declared")[0]
                             print(source + ' = line: ' + sourceline + ', warn:' + a + ' ['+source_dir[-2]+']')
                          else:
                             print(source + ' = line: ' + sourceline +' ['+source_dir[-2]+']')

                      except IndexError:
                          print(source + ' = line: ' + sourceline +' ['+source_dir[-2]+']')
              else:
                  print (" ***** Can't determine source but warning exists...")
              if debug: break
          else:
              source = Buffer[4].split(":")[0]
              sourceline = Buffer[4].split(":")[1]
              pattern = os.path.join(home_dir, "src") + "/*/"+source
              path = glob.glob(pattern)
              source_dir = path[0].split('/')
              if debug: print ("[DEBUG] source_dir :" + source_dir[-2])
              try:
                  if src_excluded.index(source_dir[-2]) :
                      warning_count += 1
                      print(Buffer[4].strip(), ' ['+source_dir[-2]+']')
              except ValueError:
                  pass

  logfile.close()

  # footer
  print ("**********************************************************************")
  print ("Warning count = " + str(warning_count))
  print ("**********************************************************************")
  return warning_count


# ---------------------------------------------------------------------------
if __name__ == "__main__":

  warno = sys.argv[1]
  try:
    home_dir = os.path.abspath(sys.argv[2])
  except IndexError:
    home_dir = ""

  sys.exit(main(warno, home_dir=home_dir))
