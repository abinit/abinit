#!/usr/bin/env python
"Check documentation and input variables"
from __future__ import division, print_function, absolute_import #unicode_literals, 

import sys
import os
import os.path
import glob
import re

def usage():
    print("\n Usage: docchk \n ")

def abinit_test_generator():
  def test_func(abenv):
     "Check documentation and input variables"
     top = abenv.apath_of("src")
     return main(abenv.home_dir)
  return {"test_func" : test_func}

def main(home_dir, verbose=False):

  home_dir = os.path.abspath(home_dir)

  # construct list of input keywords that appear in chkvars.F90
  chkvarsf90 = os.path.join(home_dir, "src/57_iovars/chkvars.F90")
  if (os.path.isfile(chkvarsf90)):
      varfile = open(chkvarsf90)
  else:
      print(" \n File ", chkvarsf90," not found! ")
      sys.exit(2)

  in_block = False
  words = []
  for line in varfile:
      if line.find("admitted variable names") > 0:
          in_block = True
      if line.find("Extra token") > 0:
          in_block = False
      if in_block == True and line.find("list_var") > 0:
          line_words=(line.split("'")[1]).split()
          for i in range(len(line_words)):
              words.append(line_words[i])

  if not words:
      print("Found empty list of words in %s " % chkvarsf90)
      print("Perhaps someone changed the format of the file?")
      print("Please modify the code in " + __file__)
      sys.exit(2)

  print( " ============================================================= ")
  print( " ABINIT Input variables: Regenerate html from abinit_vars.yml  ")
  print( " ============================================================= ")
  pathdocdir = os.path.join(home_dir, "doc")
  cmd = "cd " + pathdocdir + " ; rm -f input_variables/generated_files/allvariables.html ; python generate_doc.py > generate_doc.log"
  os.system(cmd)
  pathlogfile = os.path.join(home_dir, "doc/generate_doc.log")
  pathpymodsdir = os.path.join(home_dir, "doc/pymods")
  cmd = "cd " + pathpymodsdir + " ; python abi_check.py > abi_check.log"
  os.system(cmd)

  with open(pathlogfile) as logfile:
    for line in logfile:
        print(line)
  pathcheckfile = os.path.join(home_dir, "doc/pymods/abi_check.log")

  with open(pathcheckfile) as checkfile:
    for line in checkfile:
        print(line)

  print( " ============================================================= ")
  print( " ABINIT Input variables: Check in documentation                ")
  print( " ============================================================= ")
  varhtml = glob.glob(os.path.join(home_dir, "doc/input_variables/generated_files/var*html"))
  varallvars = glob.glob(os.path.join(home_dir, "doc/input_variables/generated_files/allvariables.html"))
  ret_code = 0
  for iwords in range(len(words)):
      deffiles = []
      for ivarhtml in range(len(varhtml)):
          with open(varhtml[ivarhtml]) as fh: varhtmldata = fh.read()
          if words[iwords] in varhtmldata:
              deffiles.append(varhtml[ivarhtml])

      if len(deffiles) > 0:
          if verbose: print("SUCCESS: ",words[iwords]," appears in ",len(deffiles)," var*html files ")
      else:
          print("FAIL: ",words[iwords]," does not appear in any var*html files ")
          ret_code += 1

      deffiles = []
      for ivarallvars in range(len(varallvars)):
          with open(varallvars[ivarallvars]) as fh: varallvarsdata = fh.read()
          if words[iwords] in varallvarsdata:
              deffiles.append(varallvars[ivarallvars])

      if len(deffiles) > 0:
          if verbose: print("SUCCESS: ",words[iwords]," appears in ",len(deffiles)," allvariables.html file as well")
      else:
          print("FAIL: ",words[iwords]," does not appear in the central allvariables.html file ")
          ret_code += 1

  print( " ============================================================= ")
  print( " ABINIT Input variables: Check in test suite                   ")
  print( " ============================================================= ")
  for iwords in range(len(words)):
      autotest = False
      for root, dirs, files in os.walk(os.path.join(home_dir, 'tests')):
          if root.find("Input")>0:
              for ifiles in range(len(files)):
                  testfilename = os.path.join(root,files[ifiles])
                  if not testfilename.endswith(".in"):
                    #print("Ignoring", testfilename)
                    continue

                  try:
                      with open(testfilename, "rt") as fh:
                        testfileinput = fh.read()
                  except Exception as exc:
                    print("FAIL: exception while opening %s\n%s" % (testfilename, str(exc)))
                    ret_code += 1
                    continue

                  if words[iwords] in testfileinput:
                      autotest = True
                      break
          if autotest:
              break

      if autotest:
          if verbose: print("SUCCESS: ",words[iwords]," appears in automatic test suite ")
      else: 
          print("FAIL: ",words[iwords]," does not appear in automatic test suite ")
          ret_code += 1

  varfile.close()
                  
  # construct list of key words appearing in anaddb input
  invars9f90 = os.path.join(home_dir, "src/77_ddb/m_anaddb_dataset.F90")
  if (os.path.isfile(invars9f90)):
      varfile = open(invars9f90)
  else:
      print(" \n File ", invars9f90," not found! ")
      sys.exit(2)

  # Scan the source and search for the calls to intagm. Parse the arguments
  # and extract the name of the variable. The prototype of intagm is:
  #    call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'brav',tread,'INT')
  re_call = re.compile(r'\s*call\s+intagm\((.+)\)\w*', re.I)

  words = []
  for line in varfile:
      m = re_call.match(line)
      if m: 
        tokens = m.group(1).split(",")
        assert len(tokens) == 9
        words.append(tokens[-3].replace("'","").replace('"',""))

  if not words:
      print( "Found empty list of words in file %s" % invars9f90)
      print( "Perhaps someone changed the format of the file?")
      print( "Please modify the code in " + __file__)
      sys.exit(2)
  #print(words)

  print(" ============================================================= ")
  print(" ANADDB Input variables: Check in documentation                ")
  print(" ============================================================= ")
  varhtml = os.path.join(home_dir, "doc/users/generated_files/help_anaddb.html")
  for iwords in range(len(words)):
      with open(varhtml) as fh: varhtmldata = fh.read()
      if words[iwords] in varhtmldata:
          if verbose: print ("SUCCESS: ",words[iwords]," appears in ",varhtml)
      else:
          print ("FAIL: ",words[iwords]," does not appear ",varhtml)
          ret_code += 1

  print( " ============================================================= ")
  print( " ANADDB Input variables: Check in test suite                   ")
  print( " ============================================================= ")
  for iwords in range(len(words)):
      autotest = False
      for root, dirs, files in os.walk(os.path.join(home_dir, 'tests')):
          if root.find("Input")>0:
              for ifiles in range(len(files)):
                  testfilename = os.path.join(root,files[ifiles])

                  if not testfilename.endswith(".in"):
                    #print("Ignoring:", testfilename)
                    continue

                  try:
                      with open(testfilename, "rt") as fh:
                        testfileinput = fh.read()
                  except Exception as exc:
                      print("FAIL: Exception while readding %s\n%s" % (testfilename, str(exc)))
                      ret_code += 1
                      continue

                  if words[iwords] in testfileinput:
                      autotest = True
                      break
          if autotest:
              break
      if autotest:
          if verbose: 
            print("SUCCESS: ",words[iwords]," appears in automatic test suite ")
      else: 
          print("FAIL: ",words[iwords]," does not appear in automatic test suite ")
          ret_code += 1

  varfile.close()

  return ret_code

if __name__ == "__main__":

  if len(sys.argv) == 1: 
    home_dir = os.path.normpath(os.path.join(os.path.dirname(sys.argv[0]), "../.."))
  else:
    home_dir = sys.argv[1] 

  exit_status = main(home_dir, verbose=False)
  sys.exit(exit_status)
