#!/usr/bin/env python
"Looking for forbidden string(s) in ABINIT/doc abinit_vars.yml and other *.yml files"
# ======================================================================
# == Python script checking if "forbidden" string(s) -- in terms of   ==
# == "abirules" -- are present in the input variable documentation    ==
# ==                                                                  ==
# == 1- The presence of a string '\n' instead of normal carriage     ==
# ==    return usually denotes that the *.yml file                    ==    
# ==    has wrongly been edited, and that                             ==
# ==    the initial correctly formatted content of some section       ==
# ==    has now become much difficult to read ...                     ==
# ==    So, should step back, and correctly edit/modify the file !    ==
# ==                                                                  ==
# ======================================================================
from __future__ import unicode_literals, division, print_function, absolute_import

import os
import re
import sys

#Activate/Deactivate tests individually
ACTIVATE_TEST1=True   

#Forbidden statements
FORBIDDEN_LIST = [  # Note: use only lowercases
"\\n",   # Forbidden '\n' instead of carriage return
]

# ---------------------------------------------------------------------------
def abinit_test_generator():
  def test_func(abenv):
     "Looking for forbidden strings in doc/*/*.yml file(s)"
     top = abenv.apath_of("doc")
     try:
       return main(top)
     except Exception:
       import sys
       raise sys.exc_info()[1] # Reraise current exception (py2.4 compliant)
  return {"test_func" : test_func}

def main(top):
  print()
  print('---------------------------------------------------------------------')
  print(' Looking for forbidden strings in doc/*/*.yml file(s):       ')
  if ACTIVATE_TEST1:
    print('  - forbidden usage of \\n instead of carriage return                 ')
  print('---------------------------------------------------------------------')

  re_docfile = re.compile("\.yml$")

  #Initialize counters
  file_total_count=0
  stat_forbidden_string_count=0
  file_forbidden_string_count=0

  #Loop over files in doc folder
  for (root, dirs, files) in os.walk(top):
    for doc in files:
      if ( re_docfile.search(doc) ):
        file_total_count +=1
        filename=os.path.join(root,doc)
        with open(filename, "rt") as fh:
          doc_data = fh.readlines()

        #Loop over lines in the file
        lineno=0
        icount_forbidden_string=0
        for line_orig in doc_data:
          lineno += 1

          #Transform line to lower case + eliminate whitespaces
          line_lower=line_orig.lower()
          line = re.sub(" ","",line_lower)

          #Look for forbidden string
          if ACTIVATE_TEST1:
            for strg in FORBIDDEN_LIST:
              if line.find(strg) != -1:
                print('  Error: %s, line %d: found \"%s\" !' % (filename,lineno,strg))
                icount_forbidden_string +=1

        #Update counters
        stat_forbidden_string_count +=icount_forbidden_string
        if icount_forbidden_string>0: file_forbidden_string_count +=1

  #Print final message
  print( '----------->')
  print( '- There are %d yml files in the complete set.' % (file_total_count))

  exit_status = ( stat_forbidden_string_count )

  if stat_forbidden_string_count==0:
    print('- No Error or Warning !')
  else:

    if stat_forbidden_string_count>0 :
      print()
      print( '>>  %d error(s) (forbidden string(s)), appearing in %d different file(s) !' % \
            (stat_forbidden_string_count,file_forbidden_string_count))
      print()
      print( '- Please replace the forbidden \\n string by the allowed carriage return ')
      print( '-  It is likely that you have incorrectly edited the file(s), and the ')
      print( '-  best strategy might be to step back, i.e. start from a pristine file ')
      print( '-  and edit it correctly. In case of abinit_vars.yml, you are advised ')
      print( '-  to use the Abivars Java Previewer with the command java -jar Abivars.jar ')
      print( '-  See the README.txt file ... ')

  return exit_status

# ---------------------------------------------------------------------------
if __name__ == "__main__":

  if len(sys.argv) == 1: 
    top = "../../../doc"
  else:
    top = sys.argv[1] 

  exit_status = main(top) 
  sys.exit(exit_status)
