#! /usr/bin/env python 
from __future__ import print_function

import sys
import os
import yaml
import re
import argparse

################################################################################

def reformat_namelist(namelist):
  """ Change the usual bibtex Lastname1, Firstname1 and Lastname2, Firstname2 and ...
      to the Phys. Rev. B style list of names.
  """
  flag_editor=0
  newnamelist=namelist
  # Very crude treatment of the appearance of the "(Eds.)" string...
  if '(Eds.)' in namelist:
    newnamelist=newnamelist.split('(Eds.)')[0]
    flag_editor=1  
  new_name=""
  names=newnamelist.split('and')
  numnames=len(names)
  for (i,one_name) in enumerate(names):
    if ',' in one_name:
      name_initials = one_name.split(',',1)
      new_name+= name_initials[1].strip()+" "+name_initials[0].strip()
      if i==len(names)-2:
        new_name+=" and "
      elif i==len(names)-1:
        if flag_editor==1:
          new_name+=" (Eds.)"
      else:
        new_name+=", "
    else:
      new_name+=one_name
  newnamelist=new_name+","
  return newnamelist


################################################################################
###############################################################################
 
# Parsing section

################################################################################

with open('abiref.bib')  as bibtex_file:
  bibtex_str = bibtex_file.read()

bibtex_str=bibtex_str.replace('optdoi','doi')
bibtex_str=bibtex_str.replace('opturl','url')
bibtex_str=bibtex_str.replace('optURI','url')
bibtex_str=bibtex_str.replace('adsurl','url')
bibtex_str=bibtex_str.replace('href','url')

bibtex_dics=[]
bibtex_items=bibtex_str.split('@')
empty=bibtex_items.pop(0)

for item in bibtex_items:

  if 'doi.org' in item:
    print(" Error : please remove the prefix http://dx.doi.org (or similar), from the doi field.")
    print(" It should start directly with the DOI information, e.g. 10.1016 ...")
    print(" This error happended while treating item :")
    print(item)
    raise

  item_dic={}
  item=item.split('{',1)
  entrytype=item[0].strip().lower()
  if not entrytype in ["article","book","incollection","phdthesis","masterthesis"]:
    print(" Not able to treat the following entrytype:",entrytype)
    raise

  #Find the ENTRYTYPE
  item_dic={'ENTRYTYPE':entrytype}

  #Find the ID
  item=item[1].split(',',1)
  id_lower=item[0].strip().lower()
  id_upper_lower=id_lower[0].upper()+id_lower[1:]
  # Check that there is a four-digit date in the ID
  m=re.search("\d{4}",id_upper_lower,flags=0)
  if m==None:
    print("The entry %s does not have the proper format, i.e. the year should be a four-digit number" % id_upper_lower)
    raise
  else:
    result=m.group() 
    if not result[0] in ["1","2"]:
      print("The entry %s does not have the proper format, i.e. the four-digit year should start with 1 or 2" % id_upper_lower)
      raise
  item_dic['ID']=id_lower[0].upper()+id_lower[1:]

  #Take care of other fields : split the different lines to be examined.
  lines=item[1].splitlines() 
  empty=lines.pop(0)
  prev=""
  flag_parenthesis=0
  for line in lines:

    if line.strip()=="":
      continue

    # Treat the case in which the previous line had ended with a parenthesis.
    if flag_parenthesis==1 and line.strip()=='}':
      #Simulate that the previous line had ended with '},'
      flag_parenthesis=0
      newline=prev+','
    else:
      newline=prev+line.strip()
    flag_parenthesis=0
    
    len_new=len(newline)

    # Takes care of lines that might be split with a carriage return. They are identified because they do not end with '},',
    # nor with '}' followed by the next line as '}' .
    # in this case, they are simply put in 'prev', and concatenated with the next line...
    if not newline[len_new-1]==',':
      # Here treat the case of ending with '}' (the flag_parenthesis is raised), or any other ending than ',', that must be continued.
      if newline[len_new-1]=='}':
        flag_parenthesis=1
      prev=newline
      continue
    else:
      newline2=newline.strip(',').strip()
      len_new=len(newline2)
      # If the line end with a comma, but is not preceeded by a parenthesis, must be continued as well.
      if not newline2[len_new-1]=='}':
        prev=newline
        continue
      prev=''
      # Now that the correct full line has been identified, separate the key and value
      split_line=newline.split('=',1)    
      key=split_line[0].strip().lower()
      value=split_line[1].strip().strip(',')
      len_value=len(value)
      value=value[1:len_value-1]
      item_dic[key]=value

  bibtex_dics.append(item_dic)

# Constitute a dictionary based on the ID as key, with formatted reference as value
# Also, constitute the body of a txt file
lines_txt=""
reference_dic={}
for (i,ref) in enumerate(bibtex_dics):
  # Very strange behaviour of python when I try to use ref["journal"] directly ?! So, use an ugly bypass.
  # Same thing for ref["volume"]
  position=0
  flag_pages=0
  flag_eprint=0
  ENTRYTYPE=""
  ID=""
  author=""
  title=""
  booktitle=""
  journal=""
  volume=""
  pages=""
  year=""
  eprint=""
  editor=""
  publisher=""
  school=""
  address=""
  for key in ref.keys():
    if key=='ENTRYTYPE':
      ENTRYTYPE=ref.values()[position].strip()
    elif key=='ID':
      ID=ref.values()[position].strip()
    elif key=='author':
      author=ref.values()[position].strip()
    elif key=='title':
      title=ref.values()[position].strip()
    elif key=='booktitle':
      booktitle=ref.values()[position].strip()
    elif key=='journal': 
      journal=ref.values()[position].strip()
    elif key=='volume':
      volume=ref.values()[position].strip()
      flag_pages+=1
    elif key=='pages':
      pages=ref.values()[position].strip()
      flag_pages+=1
    elif key=='year':
      year=ref.values()[position].strip()
    elif key=='eprint':
      eprint=ref.values()[position].strip()
      flag_eprint=1
    elif key=='editor':
      editor=ref.values()[position].strip()
    elif key=='publisher':
      publisher=ref.values()[position].strip()
    elif key=='school':
      school=ref.values()[position].strip()
    elif key=='address':
      address=ref.values()[position].strip()
    position+=1

  #DEBUG
  print("")
  print(" flag_pages,flag_eprint,ENTRYTYPE,ID,author,title,journal,volume,pages,year,eprint:")
  print(flag_pages,flag_eprint,ENTRYTYPE,ID,author,title,journal,volume,pages,year,eprint)
  print("")
  #ENDDEBUG
  # Reformat the list of authors, starting with the initials. 
  author=reformat_namelist(author)
  editor=reformat_namelist(editor)

  #At present, treat normal journals, eprints, articles in collection, books.
  #For normal journal articles, needs both volume and pages.
  formatted=""
  if ENTRYTYPE=="article" and flag_pages==2:
    formatted=' %s "%s", %s %s, %s (%s).' % (author,title,journal,volume,pages,year)
  elif ENTRYTYPE=="article" and flag_eprint==1:
    formatted=' %s "%s", %s %s (%s).' % (author,title,journal,eprint,year)
  elif ENTRYTYPE=="incollection":
    formatted=' %s "%s", in "%s", Eds. %s (%s, %s, %s), pp. %s.' % (author,title,booktitle,editor,publisher,address,year,pages)
  elif ENTRYTYPE=="phdthesis":
    formatted=' %s "%s", PhD thesis (%s, %s, %s).' % (author,title,school,address,year)
  elif ENTRYTYPE=="masterthesis":
    formatted=' %s "%s", Master thesis (%s, %s, %s).' % (author,title,school,address,year)
  elif ENTRYTYPE=="book":
    if volume!="":
      formatted=' %s "%s", vol.%s (%s, %s, %s).' % (author,title,volume,publisher,address,year)
    else:
      formatted=' %s "%s" (%s, %s, %s).' % (author,title,publisher,address,year)

  #DOI or URL is added at the end. Note that this is optional. Also, only one is mentioned, preferentially the DOI.
  try:
    doi=ref["doi"].strip()
    if doi != "":
      formatted += ' DOI: '+doi+'.'  
  except:
    try:
      url=ref["url"].strip()
      if url != "":
        formatted += ' '+url+'.'
    except:
      pass

  reference_dic[ID]=formatted

  lines_txt+= ("[%s] %s \n") %(ID,reference_dic[ID])

# Open, write and close the txt file
file_txt = 'abiref.txt'
f_txt = open(file_txt,'w')
f_txt.write(lines_txt)
f_txt.close()
print("File abiref.txt has been written ...")

# Open, write and close the yml file
file_yml = 'abiref.yml'
f_yml = open(file_yml,'w')
for ref in bibtex_dics:
  f_yml.write('-\n')
  # Will print in the order ENTRYTYPE, ID, author, title, journal, year, volume, pages, doi, then the remaining ...
  list_key={'ENTRYTYPE':0, 'ID':1, 'author':2, 'title':3, 'journal':4, 'year':5, 'volume':6, 'pages':7, 'doi':8}
  remaining=""
  text=[]
  for i in range(0,8):
    text.append("")
  for (i,key) in enumerate(ref):
    line="  "+key+" : "+ref[key]+'\n'
    try:
      number=list_key[key]
      text[number]=line
    except:
      remaining+=line
  lines=""
  # DOI is optional
  for i in range(0,8):
    if text[i] != "":
      lines+=text[i]
  lines+=remaining
  f_yml.write(lines)
f_yml.close()
print("File abiref.yml has been written ...")
print("Work done !")
################################################################################
