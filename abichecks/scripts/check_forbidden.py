#!/usr/bin/env python
"Looking for forbidden statements in ABINIT src files"
# ======================================================================
# == Python script checking if "forbidden" statements -- in terms of  ==
# == "abirules" -- are present in ABINIT source files:                ==
# ==                                                                  ==
# == 1- Access to standard output or standard error using explicit    ==
# ==    unit numbers (0, 6, 7), not using "std_out", "std_err"        ==
# ==    or "ab_out" pre-defined keywords;                             ==
# ==       wrtout.F90 and wrtout_myproc.F90 files are ignored.        ==
# ==                                                                  ==
# == 2- Use of standard error as output unit                          ==
# ==    (identified as "not recommended").                            ==
# ==                                                                  ==
# == 3- Use of "MPI_COMM_WORLD" MPI keyword instead of one of the     ==
# ==    defined communicators: xmpi_world, mpi_enreg(:)%world_comm.   ==
# ==                                                                  ==
# == 4- Explicit "allocate" or "deallocate" statements                ==
# ==    (only macros are allowed: ABI_ALLOCATE, ABI_DEALLOCATE,       ==
# ==     ABI_DATATYPE_ALLOCATE, ABI_DATATYPE_DEALLOCATE).             ==
# ==                                                                  ==
# == 5- "call" statement not placed at the start of the line          ==
# ==    (badly interpreted by abilint).                               ==
# ==                                                                  ==
# ==  Following lines will not generate an error:                     ==
# ==   Lines beginning with: "!", "if(debug)"or  "if(ab_dbg)"         ==
# ==                                                                  ==
# ==            M. Torrent - June 2011 - rev Feb 2012 - rev Apri 2014 ==
# ======================================================================
from __future__ import unicode_literals, division, print_function, absolute_import

import os
import re
import sys

from abirules_tools import find_src_dirs

#Activate/Deactivate tests individually
ACTIVATE_TEST1=True   # Forbidden write(6
ACTIVATE_TEST2=True   # Use of stderr
ACTIVATE_TEST3=True   # MPI_COMM_WORLD
ACTIVATE_TEST4=True   # Allocate/deallocate
ACTIVATE_TEST5=False  # Call at start of line

#Forbidden write statements
WRITE_FORBIDDEN_LIST = [  # Note: use only lowercases
"print *,",
"print*,",
"write(*,",
"write(6,",
"write(06,",
"write(006,",
"wrtout(6,",
"wrtout(06,",
"wrtout(006,",
"flush(6)",
"flush(06)",
"flush(006)",
"flush_unit(6)",
"flush_unit(06)",
"flush_unit(006)",
"write(7,",
"write(07,",
"write(007,",
"wrtout(7,",
"wrtout(07,",
"wrtout(007,",
"flush(7)",
"flush(07)",
"flush(007)",
"flush_unit(7)",
"flush_unit(07)",
"flush_unit(007)",
"write(0,",
"write(00,",
"write(000,",
"wrtout(0,",
"wrtout(00,",
"wrtout(000,"
]

#Not recommended write statements
WRITE_NOTRECOMMENDED_LIST = [  # Note: use only lowercases
"write(std_err,",
"wrtout(std_err,"
]

#Ignored files when looking for forbidden write access
IGNORED_WRITE_FILES = [
"wrtout.F90",
"wrtout_myproc.F90",
"m_libpaw_tools.F90",
"defs.h" # Needed to still have a full CT-QMC standalone version
]

#Forbidden MPI_COMM_WORLD keywords
COMMWORLD_FORBIDDEN_LIST = [  # Note: use only lowercases
"mpi_comm_world"
]

#Ignored files when looking for forbidden MPI_COMM_WORLD
IGNORED_COMMWORLD_FILES = [
"m_xmpi.F90",
"m_profiling_abi.F90",
"m_libpaw_mpi.F90",
"m_malloc.F90",
"defs.h"
]

#Forbidden allocate/deallocate statements
ALLOCATE_FORBIDDEN_LIST = [  # Note: use only lowercases
"deallocate",                # "deallocate" must be put before "allocate"
"allocate"
]

#Ignored files when looking for forbidden allocate/deallocate
IGNORED_ALLOCATE_FILES = [
"abi_common.h",
"libpaw.h",
"libtetra.h",
"malloc.finc",
"defs.h"
]

#List of "call" statements
CALL_STATEMENT_LIST = [  # Note: use only lowercases
"call"
]

#Ignored files when looking for forbidden call
IGNORED_CALL_FILES = [  # Note: use only lowercases
"abi_common.h",
"libpaw.h",
"libtetra.h",
"defs.h"
]

#Lines beginning with the following statements will not generate warnings
NO_ERROR_LIST = [  # Note: use only lowercases
"!",
"if(debug)",
"if(ab_dbg)"
]

def main():

  print()
  print('---------------------------------------------------------------------')
  print(' Looking for forbidden statements in ABINIT src files:               ')
  if ACTIVATE_TEST1 or ACTIVATE_TEST2:
    print('  - forbidden access to standard output/standard error               ')
  if ACTIVATE_TEST3:
    print('  - forbidden allocate/deallocate statements                         ')
  if ACTIVATE_TEST4:
    print('  - forbidden explicit MPI_COMM_WORLD communicator                   ')
  if ACTIVATE_TEST5:
    print('  - forbidden call statements not placed at the start of the line    ')
  print('---------------------------------------------------------------------')

  re_srcfile = re.compile("\.([Ff]|[Ff]90|finc|h)$")

  #Initialize counters
  file_total_count=0
  stat_forbidden_write_count=0
  file_forbidden_write_count=0
  stat_notrecommended_write_count=0
  file_notrecommended_write_count=0
  stat_forbidden_commworld_count=0
  file_forbidden_commworld_count=0
  stat_forbidden_allocate_count=0
  file_forbidden_allocate_count=0
  stat_forbidden_call_count=0
  file_forbidden_call_count=0

  # Loop over files in src folder
  for top in find_src_dirs():
      for root, dirs, files in os.walk(top):
        for src in files:
          if re_srcfile.search(src):
            file_total_count += 1
            filename = os.path.join(root,src)
            with open(filename, "rt") as fh:
              src_data = fh.readlines()

            #Loop over lines in the file
            lineno=0
            icount_forbidden_write=0
            icount_notrecommended_write=0
            icount_forbidden_commworld=0
            icount_forbidden_allocate=0
            icount_forbidden_call=0
            for line_orig in src_data:
              lineno += 1

              #Transform line to lower case + eliminate whitespaces
              line_lower=line_orig.lower()
              line = re.sub(" ","",line_lower)

              #Skip lines beginning with an authorized character
              ignored=0
              for strg in NO_ERROR_LIST:
                if line.find(strg) == 0: ignored=1

              if ignored == 0:
                #Look for forbidden write statements
                if ACTIVATE_TEST1 and (not src in IGNORED_WRITE_FILES):
                  for strg in WRITE_FORBIDDEN_LIST:
                    if line.find(strg) != -1:
                      print('  Error: %s, line %d: found \"%s\" !' % (filename,lineno,strg))
                      icount_forbidden_write +=1

                #Look for not recommended write statements
                if ACTIVATE_TEST2 and (not src in IGNORED_WRITE_FILES):
                  for strg in WRITE_NOTRECOMMENDED_LIST:
                    if line.find(strg) != -1 and src not in ["m_specialmsg.F90"]:
                      print('- Warning: %s, line %d: found \"%s\" !' % (filename,lineno,strg))
                      icount_notrecommended_write +=1

                #Look for forbidden MPI_COMM_WORLD statements
                if ACTIVATE_TEST3 and (not src in IGNORED_COMMWORLD_FILES):
                  for strg in COMMWORLD_FORBIDDEN_LIST:
                    if line.find(strg) != -1:
                      print('  Error: %s, line %d: found \"%s\" !' % (filename,lineno,strg))
                      icount_forbidden_commworld +=1

                #Look for forbidden allocate/deallocate statements
                if ACTIVATE_TEST4 and (not src in IGNORED_ALLOCATE_FILES):
                  ifound=0
                  for strg in ALLOCATE_FORBIDDEN_LIST:
                    ialloc=line.find(strg)
                    if ifound==0 and ialloc != -1:
                      ifound=1
                      if not '_'+strg in line:
                        if ialloc+len(strg)<len(line):
                          if line[ialloc-4:ialloc] != "abi_" and line[ialloc+len(strg)] == "(":
                            print('  Error: %s, line %d: found \"%s\" !' % (filename,lineno,strg))
                            icount_forbidden_allocate +=1

                #Look for forbidden call statements (not placed at the start of the line)
                if ACTIVATE_TEST5 and (not src in IGNORED_CALL_FILES):
                  for strg in CALL_STATEMENT_LIST:
                    icall=line_lower.find(" "+strg+" ")
                    if icall > 0:
                      line_before=re.sub(" ","",line_lower[:icall])
                      line_after=re.sub(" ","",line_lower[icall+3:len(line_lower)-1])
                      nothing_before=(len(line_before) == 0)
                      comment_before=(line_before.find("!") != -1)
                      label_before=(not (re.match('[0-9]+',line_before)==None))
                      quote_before=((line_before.find("'") != -1) or (line_before.find('"') != -1)) 
                      problem_before=((not nothing_before) and (not label_before) and (not quote_before))
                      problem_after=(re.match('[^\;\:\>\<\=\+\-\*\/\!\,]*(\(.*[\)\&]+|)',line_after)==None)
                      if (not comment_before) and (problem_before or problem_after):
                        print('  Error: %s, line %d: found \"%s\" not placed at the start of the line !' \
                              % (filename,lineno,strg))
                        icount_forbidden_call +=1

            #Update counters
            stat_forbidden_write_count +=icount_forbidden_write
            stat_notrecommended_write_count +=icount_notrecommended_write
            stat_forbidden_commworld_count +=icount_forbidden_commworld
            stat_forbidden_allocate_count +=icount_forbidden_allocate
            stat_forbidden_call_count +=icount_forbidden_call
            if icount_forbidden_write>0: file_forbidden_write_count +=1
            if icount_notrecommended_write>0: file_notrecommended_write_count +=1
            if icount_forbidden_commworld>0: file_forbidden_commworld_count +=1
            if icount_forbidden_allocate>0: file_forbidden_allocate_count +=1
            if icount_forbidden_call>0: file_forbidden_call_count +=1

  # Print final message
  print( '----------->')
  print( '- There are %d F90 or header files in the complete set.' % file_total_count)
  assert file_total_count

  exit_status = (stat_forbidden_write_count + stat_notrecommended_write_count + 
                 stat_forbidden_commworld_count + stat_forbidden_allocate_count +
                 stat_forbidden_call_count)

  if stat_forbidden_write_count==0 and stat_notrecommended_write_count==0 and \
     stat_forbidden_commworld_count==0 and stat_forbidden_allocate_count==0 and \
     stat_forbidden_call_count==0:
    print('- No Error or Warning !')
  else:

    if stat_forbidden_write_count > 0:
      print()
      print( '>>>  %d error(s) (forbidden write statement(s)), appearing in %d different file(s)!' % \
            (stat_forbidden_write_count,file_forbidden_write_count))
      print()
      print( '  Replace the forbidden statement(s) by allowed ones:\n')
      print( '      "write(std_out,", "write(std_err," or "write(ab_out,".\n')
      print( '  Note that std_out redirects to the ABINIT log file')
      print( '  while ab_out redirects to the ABINIT output file')
      print( '  while std_err redirects to number 0 unit.')

    if stat_notrecommended_write_count>0:
      print()
      print( '>>> %d WARNINGS (not recommended write statement(s)), appearing in %d different file(s) !' % \
            (stat_notrecommended_write_count,file_notrecommended_write_count))
      print()
      print( '  Writing to the standard error is allowed in the following cases:\n')
      print( '      - Within debugging sections of ABINIT, not activated in the production version;')
      print( '      - In case an error has been detected, causing stop to ABINIT.\n')
      print( '  In other cases, writing to std_err is not recommended, and should be avoided.')

    if stat_forbidden_commworld_count>0:
      print()
      print( '>>> %d ERRORS(s) (forbidden MPI_COMM_WORLD statement(s)), appearing in %d different file(s)!' % \
            (stat_forbidden_commworld_count,file_forbidden_commworld_count))
      print()
      print( '  Replace the MPI_COMM_WORLD forbidden statement(s) by allowed ones:\n')
      print( '      - "xmpi_world" or "mpi_enreg%world_comm".')
      print( '      - MPI_COMM_WORLD is not allowed because it may be redefined in some cases\n')

    if stat_forbidden_allocate_count > 0:
      print()
      print( '>>> %d ERROR(s) (forbidden allocate/deallocate statement(s)), appearing in %d different file(s)!\n' % \
            (stat_forbidden_allocate_count,file_forbidden_allocate_count))
      print()
      print( '   Replace the forbidden allocate statement(s) by the ABI_ALLOCATE macro defined in abi_common.h')
      print( '   Replace the forbidden deallocate statement(s) by the ABI_DEALLOCATE macro defined in abi_common.h')
      print( '   You need to add `#include "abi_common.h"` and `use m_errors` to have access to these macros.')        
      print( '   Note that the syntax of the ABI_ALLOCATE macro is not exactly the same as the allocate statement:\n')
      print( '       - only one array to be allocated in each ABI_ALLOCATE')
      print( '      - separate the array from the parenthesis and size with a comma')
      print( '      - example: instead of allocate(arr1(3,N), arr2(20)), you should write the allocations separately using:\n')
      print( '              ABI_ALLOCATE(arr1,(3,N))')
      print( '              ABI_ALLOCATE(arr2,(20))\n')
      print( '  Note that the syntax of the ABI_DEALLOCATE macro is not exactly the same as the deallocate statement:\n')
      print( '      - only one array to be deallocated in each ABI_DEALLOCATE')
      print( '      - example: instead of deallocate(arr1,arr2), you should write the deallocations separately using:\n')
      print( '              ABI_DEALLOCATE(arr1)')
      print( '              ABI_DEALLOCATE(arr2)\n')
      print( '  Finally, use ABI_MALLOC_SCALAR and ABI_FREE_SCALAR to allocate/free scalar entities.')
      print( '  and ABI_MOVE_ALLOC instead of Fortran move_alloc.')

    if stat_forbidden_call_count > 0:
      # MG: Why? Besides there are several macros doing this...
      print()
      print( '>> %d ERROR(s) (badly placed call statement(s)), appearing in %d different file(s) !' % \
            (stat_forbidden_call_count,file_forbidden_call_count))
      print()
      print( '  Please place the CALL statement at the start of the line.')
      print( '      Avoid: "statement 1 ; call subroutine()"')
      print( '      Avoid: "if (condition) call subroutine()"')

  return exit_status


if __name__ == "__main__":
  sys.exit(main())
