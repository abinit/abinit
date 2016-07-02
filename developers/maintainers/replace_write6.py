#!/usr/bin/env python
# ======================================================================
# == Python script replacing "forbidden" access to standard output    ==
# == or standard error in ABINIT source files by statements using     ==
# == "std_out" or "std_err" global variables.                         ==
# ==                                                                  ==
# ==  Unallowed access: access to standard output via an unit         ==
# ==                    number different from "std_out".              ==
# ==                    access to main output via an unit             ==
# ==                    number different from "ab_out".               ==
# ==                    access to standard output via an unit         ==
# ==                    number different from "std_err".              ==
# ==                                                                  ==
# ==  Usage: simply launch the script from the ABINIT root folder     ==
# ==                                                                  ==
# ==                                           M. Torrent - June 2011 ==
# ======================================================================
from __future__ import print_function

import os
import re

#The following statements (access to standard output)...
write6_wrong = [  # Note: use only lowercases
"print*,",
"write(*,",
"write(6,",
"write(06,",
"write(006,",
"wrtout(6,",
"wrtout(06,",
"wrtout(006,"
]
#... have to be replaced by these ones:
write6_correct = [
"write(std_out,*) ",
"write(std_out,",
"write(std_out,",
"write(std_out,",
"write(std_out,",
"wrtout(std_out,",
"wrtout(std_out,",
"wrtout(std_out,"
]

#The following statements (access to main output)...
write7_wrong = [  # Note: use only lowercases
"write(7,",
"write(07,",
"write(007,",
"wrtout(7,",
"wrtout(07,",
"wrtout(007,"
]
#... have to be replaced by these ones:
write7_correct = [
"write(ab_out,",
"write(ab_out,",
"write(ab_out,",
"wrtout(ab_out,",
"wrtout(ab_out,",
"wrtout(ab_out,"
]

#The following statements (access to standard error)...
write0_wrong = [  # Note: use only lowercases
"write(0,",
"write(00,",
"write(000,",
"wrtout(0,",
"wrtout(00,",
"wrtout(000,"
]
#... have to be replaced by these ones:
write0_correct = [
"write(std_err,",
"write(std_err,",
"write(std_err,",
"wrtout(std_err,",
"wrtout(std_err,",
"wrtout(std_err,"
]

#Pattern to be search for each member of write6_wrong
#(case insensitive, eventual whitespaces between characters)
pattern_write6_wrong=[]
for i in range(len(write6_wrong)):
  strg=re.sub(r'(?<!^)x*','\s*',write6_wrong[i])
  strg=re.sub(r'\*\*','*\*',strg)
  strg=re.sub(r'\(','\\\(',strg)
  strg='(?i)'+strg
  re.compile(strg)
  pattern_write6_wrong.append(strg)

#Pattern to be search for each member of write7_wrong
#(case insensitive, eventual whitespaces between characters)
pattern_write7_wrong=[]
for i in range(len(write7_wrong)):
  strg=re.sub(r'(?<!^)x*','\s*',write7_wrong[i])
  strg=re.sub(r'\*\*','*\*',strg)
  strg=re.sub(r'\(','\\\(',strg)
  strg='(?i)'+strg
  re.compile(strg)
  pattern_write7_wrong.append(strg)

#Pattern to be search for each member of write0_wrong
#(case insensitive, eventual whitespaces between characters)
pattern_write0_wrong=[]
for i in range(len(write0_wrong)):
  strg=re.sub(r'(?<!^)x*','\s*',write0_wrong[i])
  strg=re.sub(r'\*\*','*\*',strg)
  strg=re.sub(r'\(','\\\(',strg)
  strg='(?i)'+strg
  re.compile(strg)
  pattern_write0_wrong.append(strg)

re_srcfile = re.compile("\.([Ff]|[Ff]90|h)$")

print( )
print( '---------------------------------------------------')
print( 'Replacing forbidden access to standard output/error')
print( 'in ABINIT src files')
print( '---------------------------------------------------')

#Initialize counters
file_total_count=0
file_write_wrong_count=0

#Loop over files in src folder
#-----------------------------
for (root, dirs, files) in os.walk("src"):
  for src in files:
    if ( re_srcfile.search(src) ):
      file_total_count +=1
      filename=os.path.join(root,src)
      with open(filename, "r") as fh: 
        src_data = fh.readlines()

#     1- look for wrong write6/7/0 statements in the file
#     ---------------------------------------------------
#     Loop over lines in the file
      found_write6_wrong=False
      found_write7_wrong=False
      found_write0_wrong=False
      for line in src_data:
#       Transform line to lower case + eliminate whitespaces
        line=line.lower()
        line = re.sub(" ","",line)
#       Look for wrong statements
        for strg in write6_wrong:
          if line.find(strg) != -1:
            found_write6_wrong=True
            break
        for strg in write7_wrong:
          if line.find(strg) != -1:
            found_write7_wrong=True
            break
        for strg in write0_wrong:
          if line.find(strg) != -1:
            found_write0_wrong=True
            break
        if found_write6_wrong: break
        if found_write7_wrong: break
        if found_write0_wrong: break

#     2- if wrong statements found, replace them
#     ------------------------------------------
      if found_write6_wrong or found_write7_wrong or found_write0_wrong:
        ErrorEncountered=False
#       Write message
        print ('File %s: found wrong statements !' % (filename))
#       Open a temporary file for writing
        filenametmp=filename+'.tmp'
        try:
          filout=open(filenametmp,'w')
        except:
          print ('File %s, error: could not open tmp file !' % (filename))
          ErrorEncountered=True
        if not ErrorEncountered:
#         Loop over lines in the file
          for line in src_data:
            newline=line
#           Replace wrong statements
            for i in range(len(write6_wrong)):
              newline=re.sub(pattern_write6_wrong[i],write6_correct[i],newline)
            for i in range(len(write7_wrong)):
              newline=re.sub(pattern_write7_wrong[i],write7_correct[i],newline)
            for i in range(len(write0_wrong)):
              newline=re.sub(pattern_write0_wrong[i],write0_correct[i],newline)
#           Write the modified line
            try:
              filout.write(newline)
            except:
              print ('File %s, error: could not write into tmp file !' % (filename))
              ErrorEncountered=True
              break
#         Close the temporary file
          filout.close()
#       Replace current file by temporary file
        if not ErrorEncountered:
          try:
            os.system('mv -f '+filenametmp+' '+filename)
            file_write_wrong_count +=1
          except:
            print ('File %s, error: could not move tmp file !' % (filename))

#Print final message
#-------------------
if file_write_wrong_count==0:
  print ('No file modified !')
else:
  print ('---------------------------------------------------')
  print ('%d/%d file(s) modified !' %  (file_write_wrong_count,file_total_count))
