#! /usr/bin/env python 
from __future__ import print_function

import sys
import os
import yaml
import re
import argparse

################################################################################
###############################################################################
 
# Parsing section

################################################################################

with open('abiref.bib')  as bibtex_file:
  bibtex_str = bibtex_file.read()

bibtex_str=bibtex_str.replace('optdoi','doi')
bibtex_str=bibtex_str.replace('opturl','url')

bibtex_dics=[]
bibtex_items=bibtex_str.split('@')
empty=bibtex_items.pop(0)

for item in bibtex_items:
  item_dic={}
  item=item.split('{',1)
  item_dic={'ENTRYTYPE':item[0].strip()}
  item=item[1].split(',',1)
  item_dic['ID']=item[0].strip()
  lines=item[1].splitlines() 
  empty=lines.pop(0)
  prev=""
  for line in lines:
    newline=prev+line.strip()
    if not newline in ['}','']:
      len_new=len(newline)
      # Takes care of lines that might be split with a carriage return. They are identified because they do not end with '},' .
      # in this case, they are simply put in 'prev', and concatenated with the next line...
      if not newline[len_new-1]==',':
        prev=newline
        continue
      else:
        newline2=newline.strip(',').strip()
        len_new=len(newline2)
        if not newline2[len_new-1]=='}':
          prev=newline
          continue
        prev=''
        split_line=newline.split('=',1)    
        key=split_line[0].strip()
        value=split_line[1].strip().strip(',')
        len_value=len(value)
        value=value[1:len_value-1]
        item_dic[key]=value
  bibtex_dics.append(item_dic)

# Open, write and close the yml file
file_yml = 'abiref.yml'
f_yml = open(file_yml,'w')
for item in bibtex_dics:
  f_yml.write('-\n')
  # Will print in the order ENTRYTYPE, ID, author, title, journal, year, volume, pages, doi, then the remaining ...
  list_key={'ENTRYTYPE':0, 'ID':1, 'author':2, 'title':3, 'journal':4, 'year':5, 'volume':6, 'pages':7, 'doi':8}
  remaining=""
  text=[]
  for i in range(0,8):
    text.append("")
  for (i,key) in enumerate(item):
    line="  "+key+" : "+item[key]+'\n'
    try:
      number=list_key[key]
      text[number]=line
    except:
      remaining+=line
  lines=""
  for i in range(0,8):
    if text[i] != "":
      lines+=text[i]
  lines+=remaining
  f_yml.write(lines)
f_yml.close()
print("File abiref.yml has been written ...")
print("Work done !")
################################################################################
