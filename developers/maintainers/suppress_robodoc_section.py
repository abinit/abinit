#!/usr/bin/env python
# ======================================================================
# == Python script to remove robodoc sections, e.g. PARENTS           ==
# ==  Usage: simply launch the script from the ABINIT root folder     ==
# ==  The keyword of the section can be changed,                      ==
# ==  by replacing its keystring below                                ==
# ==                                                                  ==
# ==  X.Gonze - July 19 2022 starting from     M. Torrent - June 2011 ==
# ======================================================================
from __future__ import print_function

import os
import re

#Define keystring
keystring='PARENTS'
#keystring='CHILDREN'

re_srcfile = re.compile("\.([Ff]|[Ff]90|h)$")

print( )

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

#     1- look for '!!'+keystring in the file and identify whether it is a valid keystring section
#     ------------------------------------------------------------------------------------
#     Loop over lines in the file
      found_keystring=False
      found_keystring_section=False
      for line in src_data:
#       Eliminate whitespaces, add blank at beginning
#       print ('Line: %s' % (line))
        line = re.sub(" ","",line)
#       print ('Line: %s' % (line))
#       print(line.find("!!"+keystring+"\n"))
        if line.find("!!"+keystring+"\n") != -1:
          found_keystring=True
#         print ('File %s, found keystring statement' % (filename))
        if found_keystring:
#         print (' found_keystring is True, looking for section')
          if line.find("!!") == -1: 
            found_keystring=False
#           print ('File %s, did not find keystring section' % (filename))
#         print(line.find("!!\n"))
          if line.find("!!\n") != -1: 
            found_keystring_section=True
#           print ('File %s, found keystring section' % (filename))
            break
#         print (' found_keystring is True, this line is inside section')
        if found_keystring_section: break
#     print ('File %s, examined for keystring sections, finished search' % (filename))

#     2- if keystring statements found, suppress them
#     ----------------------------------------------
      if found_keystring_section:
        ErrorEncountered=False
#       Write message
#       print ('File %s: found at least one keystring section\n' % (filename))
#       Open a temporary file for writing
        filenametmp=filename+'.tmp'

        try:
          filout=open(filenametmp,'w')
        except:
          print ('File %s, error: could not open tmp file !' % (filename))
          ErrorEncountered=True

        if not ErrorEncountered:
#         print ('File %s, succeeded to open tmp file !' % (filename))
#         Loop over lines in the file
          found_keystring=False
          write_newline=True
          for line in src_data:

#           Eliminate whitespaces, add blank at beginning
#           print ('Line: %s' % (line))
            line_wo_blks = re.sub(" ","",line)
#           print ('Line: %s' % (line_wo_blks))

            is_keystring_line=False
            if line_wo_blks.find("!!"+keystring+"\n") != -1:
              found_keystring=True
              is_keystring_line=True
              newline=line
              write_newline=False
#             print ('File %s, found keystring statement' % (filename))
            if not found_keystring:
              write_newline=True
              newline=line
#             print ('File %s, not inside keystring section' % (filename))
            if found_keystring and not is_keystring_line:
              end_keystring=False
#             print (' found_keystring is True, looking for section')
#             Abnormal end of keystring section : write stored lines
              if line_wo_blks.find("!!") == -1:  
                end_keystring=True
                found_keystring=False
                write_newline=True
                newline=newline+line
#               print ('File %s, abnormal end of keystring section, will write stored newline' % (filename))
#             Normal end of keystring section
              if line_wo_blks.find("!!\n") != -1:  
                end_keystring=True
                found_keystring=False
                write_newline=False
#               print ('File %s, found keystring section, do not write' % (filename))
#             Inside keystring section : store the line
              if not end_keystring:
                write_newline=False
                newline=newline+line
#               print ('File %s, inside keystring section, do not write' % (filename))
           
#           Write the modified line
            if write_newline:
#             print ('File %s, will write newline' % (filename))
              try:
                filout.write(newline)
              except:
                print ('File %s, error: could not write into tmp file !' % (filename))
                ErrorEncountered=True
                break

#         Close the temporary file
          filout.close()

#         Replace current file by temporary file
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
