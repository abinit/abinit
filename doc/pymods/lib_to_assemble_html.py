#!/usr/bin/env python
""" Different functions to assemble the html files to be generated from the yml information. """

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

# Relative path from HTML files
js_path = "../../js_files/"
users_path = "../../users/"

# Path for yml and html files
bib_ori = "bibliography/origin_files"
bib_gen = "bibliography/generated_files"
help_ori = "users/origin_files"
help_gen = "users/generated_files"
invars_ori = "input_variables/origin_files"
invars_gen = "input_variables/generated_files"
topics_ori = "topics/origin_files"
topics_gen = "topics/generated_files"
theory_ori = "theory/origin_files"
theory_gen = "theory/generated_files"
tuto_ori = "tutorial/origin_files"
tuto_gen = "tutorial/generated_files"

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

  p = re.compile("\*\*([a-zA-Z0-9_ */<>.]*)\*\*")
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

def make_links(text,cur_key,allowed_link_seeds,backlinks,backlink):

  def replace_link(mymatch):
    key = mymatch.group()[2:-2]
    if key == cur_key:
      return "<b>"+key+"</b>"
    if key in allowed_link_seeds.keys():
      value=allowed_link_seeds[key]
      if "input_variable in " in value:
        # This is a link to an input variable
        varfile=value[18:]
        return '<a href="../../input_variables/generated_files/'+varfile+".html#"+key+'">'+key+'</a>'
      elif value=="characteristic":
        return '<a href="../../users/generated_files/help_abinit.html#'+str.replace(key.lower()," ","_")+'">'+key+'</a>'
      elif value=="special":
        return '<a href="../../input_variables/generated_files/specials.html#'+key+'">'+key+'</a>'
      elif value=="varfile":
        return '<a href="../../input_variables/generated_files/'+key+'.html">'+key+'</a>'
      elif value=="lesson":
        return '<a href="../../tutorial/generated_files/'+key+'.html">'+key+'</a>'
      elif value=="theorydoc":
        return '<a href="../../theory/generated_files/'+key+'.html">'+key+'</a>'
      elif value=="helpfile":
        return '<a href="../../users/generated_files/'+key+'.html">'+key[5:]+' help file</a>'
      elif value=="in_tests":
        return '<a href="../../'+key+'">&#126;abinit/'+key+'</a>'
      elif value=="allvariables":
        return '<a href="../../input_variables/generated_files/allvariables.html">'+key+'</a>'
      elif value=="bibID":
        result=get_year(key)
        if result != -9999 :
          backlinks[key]+=backlink+";;"
          return '[<a href="../../bibliography/generated_files/bibliography.html#'+key+'">'+key+'</a>]'

    return '<a href="#">[[FAKE LINK:'+key+']]</a>'

  p=re.compile("\\[\\[([a-zA-Z0-9_ */<>.]*)\\]\\]")
  if text is None:
    return ""
  new_text=p.sub(replace_link,text)

  return new_text

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

def read_yaml_file(ymlfile):
  """ Print reading message, read the file 'ymlfile' containing yml data, and store all such data in the returned object"""

  # Add echo for selected files
  if ymlfile== invars_ori+"/abinit_vars.yml":
    print("Read "+ymlfile+" as database input file for the input variables and their characteristics ...")
  elif ymlfile== bib_ori+"/bibhtml.yml":
    print("Read "+ymlfile+" as database input file for the list of generated files in the bibliography directory ...")
  elif ymlfile== topics_ori+"/default_topic.yml":
    print("Read "+ymlfile+" as database input file for the default components of the topic files ...")
  elif ymlfile== help_ori+"/helps.yml":
    print("Read "+ymlfile+" as database input file for the list of help files in the users directory ...")
  elif ymlfile== tuto_ori+"/lessons.yml":
    print("Read "+ymlfile+" as database input file for the list of lessons ...")
  elif ymlfile== topics_ori+"/list_of_topics.yml":
    print("Read "+ymlfile+" as database input file for the list of topics ...")
  elif ymlfile== topics_ori+"/tests_dirs.yml":
    print("Read "+ymlfile+" as database input file for the list of directories in which automatic test input files are present ...")
  elif ymlfile== theory_ori+"/theorydocs.yml":
    print("Read "+ymlfile+" as database input file for the list of theory documents ...")
  # This is a generated file, unlike the others
  elif ymlfile== topics_gen+"/topics_in_tests.yml":
    print("Generate file "+ymlfile+", to contain the list of automatic test input files relevant for each topic ...")
  elif ymlfile== invars_ori+"/varfiles.yml":
    print("Read "+ymlfile+" as database input file for the list of varfiles ...")
  elif "yml" in ymlfile[-3:]:
    if "topic_" in ymlfile or "lesson_" in ymlfile :
      pass
  else:
    raise KeyError("Please provide a message to be printed for the unknown file %s",ymlfile)

  with open(ymlfile, 'r') as f:
    ymlstructure = yaml.load(f);

  return ymlstructure

################################################################################

def bibtex2html(str_input):
  """ Convert the bibtex notations to html notations inside the string str """

  str=str_input

  #Subscripts
  for i in string.digits:
    string_old='$_'+i+'$'
    string_new="<sub>"+i+"</sub>"
    str=str.replace(string_old,string_new)

  #Greek letters
  list_signs=["alpha","beta","gamma","epsilon","delta","zeta","eta","theta","iota","kappa","lambda","mu","nu","xi","omicron","pi","rho","sigma","tau","upsilon","phi","chi","psi","omega"]
  list_signs_uplower=[]
  for (i,item) in enumerate(list_signs):
    list_signs_uplower.append(item[0].upper()+item[1:])
  list_signs.extend(list_signs_uplower)
  for i in list_signs:
    string_old='$\\'+i+'$'
    string_new='&'+i+';'
    str=str.replace(string_old,string_new)

  #Accented characters
  list_vowels=["a","e","i","o","u","y","A","E","I","O","U","Y"]
  list_signs_in=['"',"'","`","^","~"]
  list_signs_out=["uml","acute","grave","circ","tilde"]
  for vowel in list_vowels:
    for (i,item) in enumerate(list_signs_in):
      string_1= "{\\" + item + vowel + "}"
      string_2= "\\" + item + vowel
      string_3= "\\" + item + "{" + vowel + "}"
      #Note that the ampersand (like the backslash) is a special character in Python, so the first string is prepended with 'r' to avoid special treatment.
      string_final= r"&%s%s;" %(vowel,list_signs_out[i])
      str=str.replace(string_1,string_final)
      str=str.replace(string_2,string_final)
      str=str.replace(string_3,string_final)

  str=str.replace(r"{\~n}","&ntilde;")
  str=str.replace(r"\~n","&ntilde;")
  str=str.replace(r"\~{n}","&ntilde;")
  str=str.replace(r"{\'n}","&#324;")
  str=str.replace(r"\'n","&#324;")
  str=str.replace(r"\'{n}","&#324;")

  str=str.replace(r"{\c c}","&ccedil;")
  str=str.replace(r"\c c","&ccedil;")
  str=str.replace(r"\c{c}","&ccedil;")

  #Get rid of uneeded parentheses. One is however left with the {XYZ..} case, that should be handled with a regular expression. (TO BE DONE)
  for i in string.letters:
    string_old='{'+i+'}'
    string_new=i
    str=str.replace(string_old,string_new)
  #Here, do it on a case-by-case basis. Very unsatisfactory...
  list_signs=["ABINIT","AIP","ATOMPAW","CPU","DFT","DMFT","ELPA","ESPRESSO","GGA","GPU","GW","III","LDA","MO","PA","PAW","QE","QMR","QUANTUM","RPA","SIAM","VESTA","XML"]
  for i in list_signs:
    string_old='{'+i+'}'
    string_new=i
    str=str.replace(string_old,string_new)

  str=str.replace("--","&ndash;")

  #Suppose all remaining parenthesis are present to avoid BibTex to switch automatically from uppercase to lowercase,
  #which will not happen in HTML...
  str=str.replace('"{','"')
  str=str.replace('}"','"')

  return str

################################################################################

def assemble_html(origin_yml_files,suppl_components,dir_root,name_root,allowed_link_seeds,backlinks):
  """ Use the list of dictionaries "origin_yml_files" as well as the
      supplementary components "suppl_components", to produce html files,
      situated in dir_root directories dir_root+"/generated_files".
      The root name of the files is "name_root".
      The complementary yml information (intro, body) for each file
      is situated in dir_root+"/origin_files".
      The dictionary allowed_link_seeds allows one to set up the links to relevant keywords.
      The backlinks are accumulated, to be mentioned in the bibliography.html file.
  """

  # Store the default informations
  for i, origin_yml in enumerate(origin_yml_files):
    if origin_yml.name.strip()=="default":
      origin_yml_default=origin_yml

  # Generate each "normal" file : assemble the content, apply global transformations, then write.
  list_names=[]
  list_titles=[]
  for i, origin_yml in enumerate(origin_yml_files):
    name = origin_yml.name
    if name=="default":
      continue
    list_names.append(name)
    list_titles.append(origin_yml.keyword+' - '+origin_yml.subtitle)
  
    if name_root != "":
      full_name=name_root+"_"+name
    else:
      full_name=name

    # Try to complete the information by reading a yml file
    path_yml="%s/origin_files/%s.yml" %(dir_root,full_name)
    doc_yml={}
    if os.path.isfile(path_yml):
      print("Read "+path_yml+" to build "+full_name+".html ... ",end="")
      doc_yml=read_yaml_file(path_yml)
 
    # Try to complete the information from suppl_components
    suppl={}
    if name in suppl_components.keys():
      suppl=suppl_components[name]

    #Write a first version of the html file, in the order "header" ... up to the "end"
    #Take the info from the components "default" if there is no information on the specific document provided in the yml file.
    #Also, format specifically selected components.

    doc_html=""
    for j in ["header","title","subtitle","purpose","advice","intro","copyright","links","menu",
              "tofcontent_header","toc",
              "introduction","examples","tutorials","input_variables","input_files","references",
              "content","body",
              "links","end"]:

      item=""
      # Try to get the item from different sources
      if j in doc_yml.keys() :
        item+=doc_yml[j]
      elif j in suppl.keys() :
        item+=suppl[j]
      else:
        try:
          item+=getattr(origin_yml,j)
        except:
          pass
      item=item.strip()
      if item=="" or item=="default":
        try:
          item=getattr(origin_yml_default,j)
        except:
          pass

      # If there is nothing in the component, continue
      if item=="" or item=="default" :
        continue
      
      # Possibly apply format to selected components
      if j =="subtitle":
        item = "<h2>"+item+"</h2>\n<hr>"

      # Accumulate
      doc_html += item+"\n"

    rc=finalize_html(doc_html,origin_yml,dir_root,name_root,allowed_link_seeds,backlinks) 

  list_names.sort()
  toc_all = '<a name="list"></a>'
  toc_all += '<h3><b> Alphabetical list of all files.</b></h3>'
  cur_let_all = 'A'
  toc_all += "<p>"+cur_let_all+".&nbsp;\n"

  for (ii,name) in enumerate(list_names):
    while not (name.startswith(cur_let_all.lower()) or name.startswith(cur_let_all.upper())):
      cur_let_all = chr(ord(cur_let_all)+1)
      toc_all = toc_all + "<p>"+cur_let_all+".\n"
    if name_root != "":
      full_name=name_root+"_"+name
    else:
      full_name=name
    path_html="%s/generated_files/%s.html" %(dir_root,full_name)
    toc_all = toc_all + '<br><a href="%s"/>%s</a> [%s] &nbsp;&nbsp;\n' %(path_html,full_name,list_titles[ii])

  all_files_html=""
  spec={'users':'help files','tutorial':'lessons of the tutorial',
        'theory':'theory documents','input_variables':'varfiles','bibliography':'generated files in the bibliography directory'}
  for j in ["header","title","subtitle","copyright","links","toc_all","links","end"]:
    if j == "toc_all":
      all_files_html += toc_all
    elif j == "subtitle":
      all_files_html += 'This document lists the %s of the ABINIT package.' %(spec[dir_root])
      all_files_html += '<script type="text/javascript" src="../generic_advice.js"> </script>'
    else:
      all_files_html += getattr(origin_yml_default,j)

  name_root_all="all_files"
  rc=finalize_html(all_files_html,origin_yml_default,dir_root,name_root_all,allowed_link_seeds,backlinks)

  return "Exit assemble_html"

################################################################################

def finalize_html(doc_html,origin_yml,dir_root,name_root,allowed_link_seeds,backlinks):
  """ Final steps of the preparation of a html file : global operations and write
      The draft text of the file to be written is contained in doc_html
      Some substitutions must be done using information in origin_yml (.name, .howto, .keyword, .authors)
      The links are build using allowed_link_seeds.
      The backlinks are given back.
  """

  if name_root == "all_files":
    full_name=name_root
  elif name_root != "":
    full_name=name_root+"_"+origin_yml.name
  else:
    full_name=origin_yml.name
  backlink=' &nbsp; <a href="../../%s/generated_files/%s.html">%s</a> &nbsp; ' %(dir_root,full_name,full_name)

  doc_html=doc_html.replace("__JS_PATH__",js_path)

  if origin_yml.howto != "":
    doc_html=doc_html.replace("__HOWTO__",origin_yml.howto)

  if origin_yml.keyword != "":
    doc_html=doc_html.replace("__KEYWORD__",origin_yml.keyword)

  if origin_yml.authors != "":
    doc_html=doc_html.replace("__AUTHORS__",origin_yml.authors)

  doc_html = doku2html(make_links(doc_html,None,allowed_link_seeds,backlinks,backlink))

  #Write the finalized html file.
  path_html = "%s/generated_files/%s.html" %(dir_root,full_name)
  file_html = open(path_html,'w')
  file_html.write(doc_html)
  file_html.close()
  print("File %s written ..."%path_html)

  return "Exit finalize_html"

################################################################################
