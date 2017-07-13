#! /usr/bin/env python 
"""This script generates part of the ABINIT documentation in html for the Web site"""

from __future__ import print_function

import sys
import os
import yaml
import re
import string
import argparse

# We don't install with setup.py hence we have to add the directory [...]/abinit/doc to $PYTHONPATH
# See similar procedure in tests/runtests.py
pack_dir, x = os.path.split(os.path.abspath(__file__))
pack_dir, x = os.path.split(pack_dir)
sys.path.insert(0,pack_dir)

import doc

from doc.pymods.variables import *

debug = 0

# Relative path from HTML files
js_path = "../../js_files/"
users_path = "../../users/"

# Path for yml and html files
bib_yml = "bibliography/yml_files"
bib_gen = "bibliography/files_generated"
invars_yml = "input_variables/yml_files"
invars_html_gen = "input_variables/html_generated"
topics_yml = "topics/yml_files"
topics_work = "topics/workspace"
topics_html_gen = "topics/html_generated"

################################################################################
###############################################################################
 
# Function definitions

###############################################################################

def format_dimensions(dimensions):

  if dimensions is None:
    s = ''
  elif dimensions == "scalar":
    s = ''
  else:
    #s = str(dimensions)
    if isinstance(dimensions,list):
      s = '('
      for dim in dimensions:
        s += str(dim) + ','

      s = s[:-1]
      s += ')'
    else:
      s = str(dimensions)
 
  return s

################################################################################

def doku2html(text):

  def replace_link(mymatch):
    abivarname = mymatch.group()[2:-2]
    return "<b>"+abivarname+"</b>"

  p = re.compile("\*\*([a-zA-Z0-9_ */<>]*)\*\*")
  text2 = p.sub(replace_link,text)

  return text2

################################################################################

def format_default(defaultval):

  if defaultval is None:
    s = 'No default'
  else:
    s = "Default is "+str(defaultval)
 
  return s

################################################################################

def make_links(text,cur_abivarname,variables,characteristics,specials,backlinks,backlink):

  def replace_link(mymatch):
    abivarname = mymatch.group()[2:-2]
    if abivarname == cur_abivarname:
      return "<b>"+cur_abivarname+"</b>"
    elif abivarname in variables.keys():
      section = variables[abivarname]
      return '<a href="../../input_variables/html_generated/'+section+".html#"+abivarname+"\">"+abivarname+"</a>"
    elif abivarname in characteristics:
      return '<a href="'+users_path+'abinit_help.html#'+str.replace(abivarname.lower()," ","_")+"\">"+abivarname+"</a>"
    elif abivarname in specials:
      return '<a href="specials.html#'+abivarname+'">'+abivarname+'</a>'
    else:
      result=get_year(abivarname)
      if result != -9999 :
        backlinks[abivarname]+=backlink+";"
        return '[<a href="../../bibliography/files_generated/bibliography.html#'+abivarname+'">'+abivarname+'</a>]'
      else:
        return '<a href="#">[[FAKE LINK:'+abivarname+']]</a>'
    return mymatch.group()

  p=re.compile("\\[\\[([a-zA-Z0-9_ */<>]*)\\]\\]")
  if text is None:
    return ""
  new_text=p.sub(replace_link,text)

  return new_text

################################################################################

def read_yaml_file(ymlfile):
  """ Read the file 'ymlfile', containing yml data, and store all such data in the returned object"""

  if ymlfile== invars_yml+"/abinit_vars.yml":
    print("Use "+ymlfile+" as database input file for the input variables and their characteristics ...")
  elif ymlfile== topics_yml+"list_of_topics.yml":
    print("Use "+ymlfile+" as database input file for the list of topics ...")
  elif ymlfile== invars_yml+"/sections.yml":
    print("Use "+ymlfile+" as database input file for the list of sections ...")
  elif ymlfile== topics_yml+"/tests_dirs.yml":
    print("Use "+ymlfile+" as database input file for the list of directories in which automatic test input files are present ...")
  elif ymlfile== topics_work+"/topics_in_tests.yml":
    print("Generated file "+ymlfile+", to contain the list of automatic test input files relevant for each topic ...")

  with open(ymlfile, 'r') as f:
    ymlstructure = yaml.load(f);

  return ymlstructure

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
  names=newnamelist.split(' and ')
  numnames=len(names)
  for (i,one_name) in enumerate(names):
    if ',' in one_name:
      name_initials = one_name.split(',',1)
      new_name+= name_initials[1].strip()+" "+name_initials[0].strip()
    else:
      new_name+=one_name.strip()
    if i==len(names)-2:
      new_name+=" and "
    elif i==len(names)-1:
      if flag_editor==1:
        new_name+=" (Eds.)"
    else:
      new_name+=", "
  newnamelist=new_name+","
  return newnamelist

################################################################################

def get_year(name):
  """ Find a substring with four digits, beginning with 1 or 2, in the "name" string, and return it is succesfull. 
      Otherwise, return -9999 if there is no such string """
  m=re.search("\d{4}",name,flags=0)
  if m==None:
    result=-9999
  else:
    result=m.group()
    if not result[0] in ["1","2"]:
      result=-9999
  return result

################################################################################
################################################################################

# Parsing section, also preliminary treatment of list of topics in tests and bibliography 

################################################################################
 
variables=read_yaml_file(invars_yml+"/abinit_vars.yml")
list_of_topics=read_yaml_file(topics_yml+"/list_of_topics.yml")
sections=read_yaml_file(invars_yml+"/sections.yml")
tests_dirs=read_yaml_file(topics_yml+"/tests_dirs.yml")
bibhtml=read_yaml_file(bib_yml+"/bibhtml.yml")

with open('bibliography/abiref.bib')  as bibtex_file:
  bibtex_str = bibtex_file.read()

################################################################################
# Parse the ABINIT input files, in order to find the possible topics to which they are linked -> topics_in_tests

try :
  rm_cmd = "rm "+topics_work+"/topics_in_tests.yml"
  retcode = os.system(rm_cmd)
except :
  if debug==1 :
    print("rm "+topics_work+"/topics_in_tests.yml failed")
    print("the file was likely non existent")

for tests_dir in tests_dirs :
  grep_cmd = "grep topics input_variables/tests/%s/Input/*.in > %s/topics_in_tests.txt"%(tests_dir,topics_work)
  retcode = os.system(grep_cmd)
  if retcode == 0 :
    sed_cmd = "sed -e 's/^/- /' %s/topics_in_tests.txt >> %s/topics_in_tests.yml"%(topics_work,topics_work)
    retcode = os.system(sed_cmd)

topics_in_tests=read_yaml_file(topics_work+"/topics_in_tests.yml")
if debug==1 :
  print(" topics_in_tests :")
  print(topics_in_tests)

################################################################################
# Constitutes a list of bib items, each being a dictionary

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
  result=get_year(id_upper_lower)
  if result== -9999 :
    print("The entry %s does not have the proper format, i.e. the year should be a four-digit number starting with 1 or 2" % id_upper_lower)
    raise
  item_dic['ID']=id_upper_lower

  #Store the remaining, for later possible reordering, without any of the later treatments.
  item_dic['body']=item[1]

  item[1]=item[1].replace('optdoi','doi')
  item[1]=item[1].replace('opturl','url')
  item[1]=item[1].replace('optURI','url')
  item[1]=item[1].replace('adsurl','url')
  item[1]=item[1].replace('href','url')

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
      newline=prev+" "+line.strip()
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

################################################################################
################################################################################

# Treat first the bibliography
# Constitute a dictionary of formatted references.
# The ID is the key and the formatted reference is the value

reference_dic={}
for (i,ref) in enumerate(bibtex_dics):
  # Very strange behaviour of python when I try to use ref["journal"] directly ?! So, use an ugly bypass.
  # Same thing for ref["volume"]. In any case, one should define a class ...
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
      formatted += (' <br> DOI: <a href="https://doi.org/%s">%s</a>.') %(doi,doi)
  except:
    try:
      url=ref["url"].strip()
      if url != "":
        formatted += (' <br> URL: <a href="%s"> %s</a>.') %(url,url)
    except:
      pass

  reference_dic[ID]=formatted

################################################################################
# Order the references

bibtex_dics=sorted( bibtex_dics, key= lambda item: item['ID'])

################################################################################
# Initialize a dictionary of "back"links

backlinks=dict()
for ref in bibtex_dics:
  ID=ref["ID"]
  backlinks[ID]=" "

################################################################################
# Write a txt file, for checking purposes

lines_txt=""
for ref in bibtex_dics:
  ID=ref["ID"]
  lines_txt+= ("[%s] %s \n") %(ID,reference_dic[ID])

# Open, write and close the txt file
file_txt = bib_gen+'/abiref.txt'
f_txt = open(file_txt,'w')
f_txt.write(lines_txt)
f_txt.close()
print("File %s has been written ..." %file_txt)

################################################################################
# Write a yml file, for checking purposes

lines_yml=""
for ref in bibtex_dics:
  lines_yml+='-\n'
  # Will print in the order ENTRYTYPE, ID, author, title, journal, year, volume, pages, doi, then the remaining (excluding the "body") ...
  list_key={'ENTRYTYPE':0, 'ID':1, 'author':2, 'title':3, 'journal':4, 'year':5, 'volume':6, 'pages':7, 'doi':8}
  remaining=""
  text=[]
  for i in range(0,8):
    text.append("")
  for (i,key) in enumerate(ref):
    if key=="body":
      continue
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
  lines_yml+=lines

# Open, write and close the yml file
file_yml = bib_gen+'/abiref.yml'
f_yml = open(file_yml,'w')
f_yml.write(lines_yml)
f_yml.close()
print("File %s has been written ..." %file_txt)

################################################################################
################################################################################

# Generate the different "section" files (var*.html, specials.html, allvariables.html)

################################################################################
# Initialization to constitute the body of the specials and var*.html files

all_vars = dict()
all_contents = dict()
for i, section_info in enumerate(sections):
  section = section_info.name
  all_vars[section] = []
  all_contents[section]= "<br><br><br><br><hr>\n"

# Create a dictionary that gives the section for each abivarname
list_all_vars = dict()
for var in variables:
  list_all_vars[var.abivarname] = var.section

################################################################################
# Constitute the body of information for the special parameters, stored for the appropriate section in all_contents[section]

cur_specials = []
for (specialkey,specialval) in list_specials:
  cur_specials.append(specialkey)

for (speckey, specval) in list_specials:
  backlink= ' <a href="../../%s/specials.html#%s">%s</a> ' %(invars_html_gen,speckey,speckey)
  cur_content = "<br><font id=\"title\"><a name=\""+speckey+"\">"+speckey+"</a></font>\n"
  cur_content += "<br><font id=\"text\">\n"
  cur_content += "<p>\n"+doku2html(make_links(specval,speckey,list_all_vars,list_chars,cur_specials,backlinks,backlink))+"\n"
  cur_content += "</font>"
  cur_content += "<br><br><a href=#top>Go to the top</a>\n"
  cur_content += "<B> | </B><a href=\"allvariables.html#top\">Complete list of input variables</a><hr>\n"
  #
  all_contents["specials"] = all_contents["specials"] + cur_content + "\n\n"

################################################################################
# Constitute the body of information for all variables, stored for the appropriate section in all_contents[section]

for i, var in enumerate(variables):
  # Constitute the body of information related to one input variable
  section = var.section
  all_vars[section].append([var.abivarname,var.mnemonics])
  cur_content = ""
  backlink=' <a href="../../%s/%s.html#%s">%s</a> ' %(invars_html_gen,section,var.abivarname,var.abivarname)

  try:
    # Title
    cur_content += "<br><font id=\"title\"><a name=\""+var.abivarname+"\">"+var.abivarname+"</a></font>\n"
    # Mnemonics
    cur_content += "<br><font id=\"mnemonics\">Mnemonics: "+var.mnemonics+"</font>\n"
    # Characteristics
    if var.characteristics is not None:
      chars = ""
      for chs in var.characteristics:
        chars += chs+", "
      chars = chars[:-2]
      cur_content += "<br><font id=\"characteristic\">Characteristic: "+make_links(chars,var.abivarname,list_all_vars,list_chars,cur_specials,backlinks,backlink)+"</font>\n"
    else:
      cur_content += "<br><font id=\"characteristic\">Characteristic: </font>\n"
    # Topics
    try:
      if var.topics is not None :
        cur_content += "<br><font id=\"characteristic\">Mentioned in \"How to\": "
        vartopics=var.topics
        topics_name_class = vartopics.split(',')
        for i, topic_name_class in enumerate(topics_name_class):
          name_class = topic_name_class.split('_')
          cur_content += '<a href="../../'+topics_html_gen+'/topic_'+name_class[0].strip()+'.html">'+name_class[0].strip()+'</a> '
        cur_content += "</font>\n"
    except:
      if debug==1 :
        print(" No topic_class for abivarname "+var.abivarname)
    # Variable type, including dimensions
    cur_content += "<br><font id=\"vartype\">Variable type: "+var.vartype
    if var.dimensions is not None:
      cur_content += make_links(format_dimensions(var.dimensions),var.abivarname,list_all_vars,list_chars,cur_specials,backlinks,backlink)
    if var.commentdims is not None and var.commentdims != "":
      cur_content += " (Comment: "+make_links(var.commentdims,var.abivarname,list_all_vars,list_chars,cur_specials,backlinks,backlink)+")"
    cur_content += "</font>\n" 
    # Default
    cur_content += "<br><font id=\"default\">"+make_links(format_default(var.defaultval),var.abivarname,list_all_vars,list_chars,cur_specials,backlinks,backlink)
    if var.commentdefault is not None and var.commentdefault != "":
      cur_content += " (Comment: "+make_links(var.commentdefault,var.abivarname,list_all_vars,list_chars,cur_specials,backlinks,backlink)+")"
    cur_content += "</font>\n" 
    # Requires
    if var.requires is not None and var.requires != "":
      cur_content += "<br><br><font id=\"requires\">\nOnly relevant if "+doku2html(make_links(var.requires,var.abivarname,list_all_vars,list_chars,cur_specials,backlinks,backlink))+"\n</font>\n"
    # Excludes
    if var.excludes is not None and var.excludes != "":
      cur_content += "<br><br><font id=\"excludes\">\nThe use of this variable forbids the use of "+doku2html(make_links(var.excludes,var.abivarname,list_all_vars,list_chars,cur_specials,backlinks,backlink))+"\n</font>\n"
    # Text
    cur_content += "<br><font id=\"text\">\n"
    cur_content += "<p>\n"+doku2html(make_links(var.text,var.abivarname,list_all_vars,list_chars,cur_specials,backlinks,backlink))+"\n"
    # End the section for one variable
    cur_content += "</font>\n\n"
    cur_content += "<br><br><a href=#top>Go to the top</a>\n"
    cur_content += "<B> | </B><a href=\"allvariables.html#top\">Complete list of input variables</a><hr>\n"
    #
    all_contents[section] = all_contents[section] + cur_content + "\n\n"
  except AttributeError as e:
    print(e)
    print('For variable : ',abivarname)

################################################################################
# Generate the files that document all the variables (all such files : var* as well as all and special).

# Store the default informations
for i, section_info in enumerate(sections):
  if section_info.name.strip()=="default":
    section_info_default=section_info

# Generate each "normal" section file : build the missing information (table of content), assemble the content, apply global transformations, then write.
for i, section_info in enumerate(sections):
  section = section_info.name
  if section=="default":
    continue

  #Generate the body of the table of content 
  cur_let = 'A'
  toc_body = " <br>"+cur_let+".\n"
  if section=="allvariables":
    for i, var in enumerate(variables):
      while not var.abivarname.startswith(cur_let.lower()):
        cur_let = chr(ord(cur_let)+1)
        toc_body += " <p>"+cur_let+".\n"
      abivarname=var.abivarname
      if var.characteristics is not None and '[[INTERNAL_ONLY]]' in var.characteristics:
        abivarname = '%'+abivarname
      curlink = " <a href=\""+var.section+".html#"+var.abivarname+"\">"+abivarname+"</a>&nbsp;&nbsp;\n"
      toc_body += curlink
  elif section == "specials":
    for (speckey, specval) in list_specials:
      while not speckey.lower().startswith(cur_let.lower()):
        cur_let = chr(ord(cur_let)+1)
        toc_body += " <br>"+cur_let+".\n"
      curlink = " <a href=\"#"+speckey+"\">"+speckey+"</a>&nbsp;&nbsp;\n"
      toc_body += curlink
  else:
    for abivarname,defi in all_vars[section]:
      while not abivarname.startswith(cur_let.lower()):
        cur_let = chr(ord(cur_let)+1)
        toc_body += " <br>"+cur_let+".\n"
      curlink = " <a href=\"#"+abivarname+"\">"+abivarname+"</a>&nbsp;&nbsp;\n"
      toc_body += curlink
  toc_body += "\n"

  #Write a first version of the html file, in the order "header" ... up to the "end"
  #Take the info from the section "default" if there is no information on the specific section provided in the yml file.
  sectionhtml=""
  for j in ["header","title","subtitle","purpose","advice","copyright","links","menu","tofcontent_header","tofcontent_body","content","links","end"]:
    if j == "tofcontent_body":
      sectionhtml += toc_body
    elif j == "content":
      sectionhtml += all_contents[section]
    else:
      extract_j=getattr(section_info,j) 
      if extract_j.strip() == "":
        sectionhtml += getattr(section_info_default,j) 
      else:
        sectionhtml += extract_j
    sectionhtml += "\n"

  #Global operations on the tentative html file.
  sectionhtml=sectionhtml.replace("__JS_PATH__",js_path)
  sectionhtml=sectionhtml.replace("__KEYWORD__",section_info.keyword)
  backlink=' <a href="../../%s/%s.html#%s">%s<\a>' %(invars_html_gen,section,var.abivarname,var.abivarname)
  sectionhtml = doku2html(make_links(sectionhtml,None,list_all_vars,list_chars,cur_specials,backlinks,backlink))

  #Write the finalized html file.
  file_cur = invars_html_gen+'/'+section+'.html'
  f_cur = open(file_cur,'w')
  f_cur.write(sectionhtml)
  f_cur.write("\n")
  f_cur.close()
  print("File %s written ..."%file_cur )

################################################################################
################################################################################

# Generate the different "topic" files

################################################################################
# Constitute the section "Related input variables" for all topic files. 
# This section in input variables is stored, for each topic_name, in topic_invars[topic_name]

topic_invars = dict()
found = dict()

for topic_name in list_of_topics:
  topic_invars[topic_name] = ""

for (tclasskey, tclassval) in list_topics_class:

  if debug == 1:
    print("\nWork on "+tclasskey+"\n")

  for topic_name in list_of_topics:
    found[topic_name] = 0

  for i, var in enumerate(variables):
    try:
      if var.topics is not None:
        topics_name_class = var.topics.split(',')
        for i, topic_name_class in enumerate(topics_name_class):
          name_class = topic_name_class.split('_')
          if tclasskey==name_class[1].strip() :
            topic_name=name_class[0].strip()
            if found[topic_name]==0 :
              topic_invars[topic_name] += "<p>"+tclassval+":<p>"
              found[topic_name] = 1
            abivarname=var.abivarname
            if var.characteristics is not None and '[[INTERNAL_ONLY]]' in var.characteristics:
              abivarname = '%'+abivarname
            topic_invars[topic_name] += '... <a href="../../'+invars_html_gen+'/'+var.section+'.html#'+var.abivarname+'">'+abivarname+'</a>   '
            topic_invars[topic_name] += "["+var.mnemonics+"]<br>\n"
    except:
      if debug==1 :
       print(" No topics for abivarname "+var.abivarname) 

################################################################################
# Constitute the section "Selected input files" for all topic files.
# This section on input files in topics  is stored, for each topic_name, in topic_infiles[topic_name]

topic_infiles = dict()

for topic_name in list_of_topics:
  topic_infiles[topic_name] = ""
 
# Create a dictionary to contain the list of tests for each topic
inputs_for_topic = dict()
for str in topics_in_tests:
  str2 = str.split(':')
  listt=str2[1]
  str_topics = listt[listt.index('=')+1:]
  list_topics = str_topics.split(',')
  for topic in list_topics:
    topic = topic.strip()
    if topic not in inputs_for_topic.keys():
      inputs_for_topic[topic] = []
    inputs_for_topic[topic].append(str2[0])

if debug==1 :
  print(inputs_for_topic)

for i, topic_name in enumerate(inputs_for_topic):
  if topic_name in list_of_topics:
    tests=inputs_for_topic[topic_name]
    # Constitute a dictionary for each directory of tests
    dir = dict()
    for test in tests:
      file_split=test.split('/')
      dirname=file_split[2] 
      testname=file_split[4] 
      if dirname not in dir.keys():
        dir[dirname] = []
      dir[dirname].append(testname)
    for dirname, testnames in dir.items():
      line="<p> tests/"+dirname+"/Input: "
      for testname in testnames:
        line+="<a href=\"../tests/"+dirname+"/Input/"+testname+"\">"+testname+"</a> \n"
      topic_infiles[topic_name]+= line
    topic_infiles[topic_name] += "<br>\n"

################################################################################
# Initiate the table of content of the file all_topics.html

toc_all = '<a name="list"></a>'
toc_all += '<h3><b> Alphabetical list of all "How to ?" documentation topics.</b></h3>'

cur_let_all = 'A'
toc_all += "<p>"+cur_let_all+".&nbsp;\n"

################################################################################
# Assemble the "topic" files 

print("Will use file yml_files/default_topic.yml as default for all topic files ... ")
default_topic_yml=read_yaml_file(topics_yml+"/default_topic.yml")
default_topic=default_topic_yml[0]

# For each "topic" file
for topic_name in list_of_topics:
  f_topic="yml_files/topic_"+topic_name+".yml"
  print("Will use file "+f_topic+" to initiate the topic "+topic_name+" ... ",end="")
  topic_yml=read_yaml_file(topics_yml+"/topic_"+topic_name+".yml")
  topic=topic_yml[0]

  #Mention it in the table of content of the file all_topics.html
  while not (topic_name.startswith(cur_let_all.lower()) or topic_name.startswith(cur_let_all.upper())):
    cur_let_all = chr(ord(cur_let_all)+1)
    toc_all = toc_all + "<p>"+cur_let_all+".\n"
  toc_all = toc_all + "<br><a href=\"topic_"+ topic_name + ".html\">" + topic_name + "</a> [How to "+topic.howto+"] &nbsp;&nbsp;\n"

  #Generate the table of content
  item_toc=0
  item_list=[]
  title={ "introduction":"Introduction." , "examples":"Example(s)", "tutorials":"Related tutorials." , "input_variables":"Related input variables." , "input_files":"Selected input files."}
  sec_number={}
  toc=" <h3><b>Table of content: </b></h3> \n <ul> "
  for j in ["introduction","examples","tutorials","input_variables","input_files"] :
    sec_number[j]="0"
    try :
      extract_j=getattr(topic,j).strip()
    except :
      extract_j=""
    if (extract_j != "" and extract_j!= "default") or (j=="input_variables" and topic_invars[topic_name]!="") or (j=="input_files" and topic_infiles[topic_name]!=""):
      item_toc += 1
      item_num="%d" % item_toc
      sec_number[j]=item_num
      toc += '<li><a href="topic_'+topic_name+'.html#'+item_num+'">'+item_num+'</a>. '+title[j]
      
  toc+= "</ul>"

  #Generate a first version of the html file, in the order "header" ... up to the "end"
  #Take the info from the section "default" if there is no information on the specific section provided in the yml file.
  topic_html=""
  for j in ["header","title","subtitle","copyright","links","toc","introduction","examples","tutorials","input_variables","input_files","links","end"]:
    if j == "toc":
      topic_html += toc
    elif j == "input_variables":
      if sec_number[j]!="0" :
        topic_html+= '\n&nbsp; \n<HR ALIGN=left> \n<a name=\"'+sec_number[j]+'\">&nbsp;</a>\n<h3><b>'+sec_number[j]+'. '+title[j]+'</b></h3>\n\n\n'
        topic_html+= topic_invars[topic_name]
    elif j == "input_files":
      if sec_number[j]!="0" :
        topic_html+= '\n&nbsp; \n<HR ALIGN=left> \n<a name=\"'+sec_number[j]+'\">&nbsp;</a>\n<h3><b>'+sec_number[j]+'. '+title[j]+'</b></h3>\n\n\n'
        topic_html+= "The user can find some related example input files in the ABINIT package in the directory /tests, or on the Web:\n"
        topic_html+= topic_infiles[topic_name]
    else:
      extract_j=getattr(topic,j).strip()
      if extract_j == "" or extract_j== "default" :
        try:
          topic_html+= getattr(default_topic,j)
        except:
          pass
      else:
        if j in title.keys():
          topic_html+= '\n&nbsp; \n<HR ALIGN=left> \n<a name=\"'+sec_number[j]+'\">&nbsp;</a>\n<h3><b>'+sec_number[j]+'. '+title[j]+'</b></h3>\n\n\n'
        topic_html += extract_j
    topic_html += "\n"

  #Global operations on the tentative html file.
  topic_html=topic_html.replace("__JS_PATH__",js_path)
  topic_html=topic_html.replace("__HOWTO__",topic.howto)
  topic_html=topic_html.replace("__KEYWORD__",topic.keyword)
  backlink=' <a href="../../%s/topic_%s.html">%s</a>' %(topics_html_gen,topic_name,topic_name)

  topic_html = doku2html(make_links(topic_html,None,list_all_vars,list_chars,cur_specials,backlinks,backlink))

  # Open, write and close the file
  file_topic = topics_html_gen+'/topic_'+topic_name+'.html'
  f_topic = open(file_topic,'w')
  f_topic.write(topic_html)
  f_topic.close()
  print("File %s written ..."%file_topic )

################################################################################
# Generate the file with the list of names of different "topic" files

#Generate a first version of the html file, in the order "header" ... up to the "end"
#Take the info from the section "default" if there is no information on the specific section provided in the yml file.
all_topics_html=""
for j in ["header","title","subtitle","copyright","links","toc_all","links","end"]:
  if j == "toc_all":
    all_topics_html += toc_all
  elif j == "subtitle":
    all_topics_html += "<h2>Complete list.</h2><hr>"
    all_topics_html += 'This document lists the names of all "How to ?" documentation topics for the abinit package.'
    all_topics_html += '<script type="text/javascript" src="../generic_advice.js"> </script>'
  else:
    all_topics_html += getattr(default_topic,j)

#Global operations on the tentative html file.
all_topics_html=all_topics_html.replace("__JS_PATH__",js_path)
all_topics_html=all_topics_html.replace("__KEYWORD__",default_topic.keyword)
backlink='<a href="../../%s/all_topics.html">all_topics.html</a>' %(topics_html_gen)
all_topics_html = doku2html(make_links(all_topics_html,None,list_all_vars,list_chars,cur_specials,backlinks,backlink))

# Open, write and close the file
file_html = topics_html_gen+'/all_topics.html'
f_html = open(file_html,'w')
f_html.write(all_topics_html)
f_html.close()
print("File %s written ..."%file_html )

################################################################################
################################################################################

# Come back to the bibliography

################################################################################
# Treat the links within the "introduction" of the acknowledgment section first.
backlink= ' <a href="../../%s/acknowledgments.html">acknowledgments.html</a> ' %(bib_gen)
for i, bibhtml_info in enumerate(bibhtml):
  if bibhtml_info.name.strip()=="acknowledgments":
    bibhtml_intro=bibhtml_info.introduction
    bibhtml_ack_intro = doku2html(make_links(bibhtml_intro,None,list_all_vars,list_chars,cur_specials,backlinks,backlink))

################################################################################
# Write an ordered bib file, that allows to update the original one.
# Constitute also the content of the ABINIT bibtex and bibliography files.

bib_content=dict()
bib_content['bibtex']=""
bib_content['bibliography']=""
bib_content['acknowledgments']=""
lines_txt=""
cur_let = 'A'
alphalinks="<hr> Goto "
for i in string.ascii_uppercase:
  alphalinks+=('<a href=#%s>%s</a> ')%(i,i)
bib_content['bibliography']+=('<a id="%s"></a>')%(cur_let)+alphalinks+('<hr><hr><h2>%s</h2>')%(cur_let)
for ref in bibtex_dics:
  entrytype=ref["ENTRYTYPE"]
  ID=ref["ID"]
  list_backlinksID=backlinks[ID].split(";")
  for (i,link) in enumerate(list_backlinksID):
     list_backlinksID[i]=link.strip()
  backlinksID=set(list_backlinksID)
  line=("@%s{%s,%s") %(entrytype,ID,ref['body'])
  lines_txt+= line
  bib_content['bibtex']+= ('<hr><a id="%s">%s</a> \n <pre>' ) %(ID,ID)
  bib_content['bibtex']+= line+'</pre> \n'
  if ID[0]>cur_let:
    cur_let=ID[0]
    bib_content['bibliography']+=('<a id="%s"></a>')%(cur_let)+alphalinks+('<hr><hr><h2>%s</h2>')%(cur_let)
  bib_content['bibliography']+= ('<hr><a id="%s">[%s]</a> (<a href="./bibtex.html#%s">bibtex</a>)\n <br> %s \n') %(ID,ID,ID,reference_dic[ID])
  nlink=0
  for link in backlinksID:
    if len(link)!=0:
      if nlink==0: 
        bib_content['bibliography']+= "<br> Referred to in " 
        nlink=1
      bib_content['bibliography']+= link

# Open, write and close the txt file
file_txt = bib_gen+'/ordered_abiref.bib'
f_txt = open(file_txt,'w')
f_txt.write(lines_txt)
f_txt.close()
print("File %s has been written ..." %file_txt)

################################################################################
#Global operation on the bibliography html file : conversion from bibtex notation to html notation.
#This should cover most of the cases. 

#Subscripts
for i in string.digits:
  string_old='$_'+i+'$'
  string_new="<sub>"+i+"</sub>"
  bib_content['bibliography']=bib_content['bibliography'].replace(string_old,string_new)

#Greek letters
list_signs=["alpha","beta","gamma","epsilon","delta","zeta","eta","theta","iota","kappa","lambda","mu","nu","xi","omicron","pi","rho","sigma","tau","upsilon","phi","chi","psi","omega"]
list_signs_uplower=[]
for (i,item) in enumerate(list_signs):
  list_signs_uplower.append(item[0].upper()+item[1:])
list_signs.extend(list_signs_uplower)
for i in list_signs:
  string_old='$\\'+i+'$'  
  string_new='&'+i+';'  
  bib_content['bibliography']=bib_content['bibliography'].replace(string_old,string_new)

#Accented characters
#Note that the backslash is a special character in Python, so the first string is prepended with 'r' to avoid special treatment.
list_vowels=["a","e","i","o","u","y","A","E","I","O","U","Y"]
list_signs_in=['"',"'","`","^","~"]
list_signs_out=["uml","acute","grave","circ","tilde"]
for vowel in list_vowels:
  for (i,item) in enumerate(list_signs_in):
    string_1= "{\\" + item + vowel + "}"
    string_2= "\\" + item + vowel
    string_3= "\\" + item + "{" + vowel + "}" 
    string_final= r"&%s%s;" %(vowel,item)
    bib_content['bibliography']=bib_content['bibliography'].replace(string_1,string_final)
    bib_content['bibliography']=bib_content['bibliography'].replace(string_2,string_final)
    bib_content['bibliography']=bib_content['bibliography'].replace(string_3,string_final)

bib_content['bibliography']=bib_content['bibliography'].replace(r"{\~n}","&ntilde;")
bib_content['bibliography']=bib_content['bibliography'].replace(r"\~n","&ntilde;")
bib_content['bibliography']=bib_content['bibliography'].replace(r"\~{n}","&ntilde;")
bib_content['bibliography']=bib_content['bibliography'].replace(r"{\'n}","&#324;")
bib_content['bibliography']=bib_content['bibliography'].replace(r"\'n","&#324;")
bib_content['bibliography']=bib_content['bibliography'].replace(r"\'{n}","&#324;")

bib_content['bibliography']=bib_content['bibliography'].replace(r"{\c c}","&ccedil;")
bib_content['bibliography']=bib_content['bibliography'].replace(r"\c c","&ccedil;")
bib_content['bibliography']=bib_content['bibliography'].replace(r"\c{c}","&ccedil;")

#Get rid of uneeded parentheses. One is however left with the {XYZ..} case, that should be handled with a regular expression. (TO BE DONE)
for i in string.letters:
  string_old='{'+i+'}'
  string_new=i
  bib_content['bibliography']=bib_content['bibliography'].replace(string_old,string_new)
#Here, do it on a case-by-case basis. Very unsatisfactory...
list_signs=["ABINIT","AIP","ATOMPAW","CPU","DFT","DMFT","ELPA","ESPRESSO","GGA","GPU","GW","III","LDA","MO","PA","PAW","QE","QMR","QUANTUM","RPA","SIAM","VESTA","XML"]
for i in list_signs:
  string_old='{'+i+'}'
  string_new=i
  bib_content['bibliography']=bib_content['bibliography'].replace(string_old,string_new)

bib_content['bibliography']=bib_content['bibliography'].replace("--","&ndash;")

#Suppose all remaining parenthesis are present to avoid BibTex to switch automatically from uppercase to lowercase,
#which will not happen in HTML...
bib_content['bibliography']=bib_content['bibliography'].replace('"{','"')
bib_content['bibliography']=bib_content['bibliography'].replace('}"','"')

################################################################################
# Generate the html files in the bibliography directory

# Store the default informations
for i, bibhtml_info in enumerate(bibhtml):
  if bibhtml_info.name.strip()=="default":
    bibhtml_info_default=bibhtml_info

# Generate each bibhtml file : build the missing information (table of content), assemble the content, apply global transformations, then write.
for i, bibhtml_info in enumerate(bibhtml):
  bibname = bibhtml_info.name
  if bibname =="default":
    continue

 #Write a first version of the html file, in the order "header" ... up to the "end"
  #Take the info from the section "default" if there is no information on the specific section provided in the yml file.
  bibhtml=""
  for j in ["header","title","subtitle","copyright","links","introduction","content","links","end"]:
    if j == "content":
      bibhtml += bib_content[bibname]
    elif j == "introduction" and bibname =="acknowledgments" :
      bibhtml += bibhtml_ack_intro
    else:
      extract_j=getattr(bibhtml_info,j)
      if extract_j.strip() == "":
        bibhtml += getattr(bibhtml_info_default,j)
      else:
        bibhtml += extract_j
    bibhtml += "\n"

  #Global operations on the tentative html file.
  bibhtml=bibhtml.replace("__JS_PATH__",js_path)
  bibhtml=bibhtml.replace("__KEYWORD__",bibhtml_info.keyword)
  backlink= ' <a href="../../%s/specials.html#%s">%s</a> ' %(bib_gen,bibname,bibname)
  bibhtml = doku2html(make_links(bibhtml,None,list_all_vars,list_chars,cur_specials,backlinks,backlink))

  #Write the finalized html file.
  file_cur = bib_gen+'/'+bibname+'.html'
  f_cur = open(file_cur,'w')
  f_cur.write(bibhtml)
  f_cur.write("\n")
  f_cur.close()
  print("File %s written ..."%file_cur )

################################################################################
#Global operation on the bibliography html file : conversion from bibtex notation to html notation.
#replace_list=[ "/{\\\'e}/&eacute;/",
#               "/\\\'e/&eacute;/",
#               "/\\\'{e}/&eacute;/"]
#os_cmd= "cp %s/bibliography.html %s/work0" %(bib_gen,bib_gen)
#retcode = os.system(os_cmd)
#for (i,item) in enumerate(replace_list):
#  sed_cmd = "sed -e '%s' %s/work%s > %s/work%s " %(item,bib_gen,i,bib_gen,i+1)
#  #DEBUG
#  print("")
#  print(" Will work on i, item:",i,item)
#  print(" sed_cmd :",sed_cmd)
#  print("")
#  #ENDDEBUG
#  retcode = os.system(sed_cmd)
#os_cmd= "cp %s/work%s %s/bibliography.html" %(bib_gen,i+1,bib_gen)
#retcode = os.system(os_cmd)
#os_cmd= "rm %s/work*" %(bib_gen)
#retcode = os.system(os_cmd)

################################################################################

print("Work done !")

################################################################################
