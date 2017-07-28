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

################################################################################
###############################################################################

# Function definitions

###############################################################################

def read_yaml(path_ymlfile):

  with open(path_ymlfile, 'r') as f:
    try:
      yml_content = yaml.load(f)
    except:
      print("\n\n... ERROR ...\n[Complement of information from generate_doc.py]")
      print("Look for ... a forbidden ':' sign in the line and file mentioned below,")
      print("      or ... an incorrect indentation in the line and file mentioned below.\n")
      raise
    return yml_content

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


def format_default(defaultval):

  if defaultval is None:
    s = 'No default'
  else:
    s = "Default is "+str(defaultval)

  return s

################################################################################

def make_links(text,cur_key,allowed_link_seeds,backlinks,backlink):
  """ Interpreter for the address contained in [[...]], following dokuwiki conventions,
      with translation to HTML links. The allowed link seeds are given.
      Exceptions to doku conventions : 
      - the cur_key matching yields simple bold emphasis of the text ;
      - bibliographical references keep one pair of square brackets.
  """

  def replace_link(mymatch):
    dokukey = mymatch.group()[2:-2].strip()
    if cur_key != None:
      if dokukey == cur_key.strip():
        return "<b>"+dokukey+"</b>"
    #Extract the four possible parts of a dokukey, with separators :, # and |
    #First, the namespace
    if ':' in dokukey:
      p1234=dokukey.split(':',1)
      namespace=p1234[0].strip()
      p234=p1234[1].strip()
    else:
      namespace=""
      p234=dokukey
    #Then, the text
    if '|' in p234:
      p234_split=p234.split('|',1) 
      p23=p234_split[0].strip()
      text=p234_split[1].strip()
    else:
      p23=p234
      text=p234
    #Finally, the key (often a filename, but not for input variables) and section
    if '#' in p23:
      p23_split=p23.split('#',1)
      key=p23_split[0].strip()
      section=p23_split[1].strip()
    else:
      key=p23
      section=""

    if namespace=="":
      linkseed=key
    else:
      linkseed=namespace+'_'+key

    #The allowed namespaces are (varfile is an exception, see later):
    dic_namespaces={"lesson":"tutorial/generated_files",
                    "theorydoc":"theory/generated_files",
                    "help":"users/generated_files",
                    "topic":"topics/generated_files",
                    "varset":"input_variables/generated_files",
                    "varfile":"input_variables/generated_files"}

    #Actually make the selection on the linkseed at present ... which allows the varfile , which is not really a namespace ...
    if linkseed in allowed_link_seeds.keys():
      value=allowed_link_seeds[linkseed]

      #Treat first the allowed namespaces
      if value in dic_namespaces.keys():
        dir=dic_namespaces[value]

        #Specific formatting treatment
        if value=="help" and namespace=="":
          text=key[5:]+' help file'

        return '<a href="../../%s/%s.html#%s">%s</a>' %(dir,linkseed,section,text)

      #Treat everything else
      elif "input_variable in " in value:
        # This is a link to an input variable
        filename=value[18:]
        return '<a href="../../input_variables/generated_files/%s.html#%s">%s</a>' %(filename,key,text)
      elif value=="characteristic":
        return '<a href="../../users/generated_files/help_abinit.html#%s">%s</a>' %(key,text)
      elif value=="in_tests":
        return '<a href="../../%s">&#126;abinit/%s</a>' %(key,key)
      elif value=="bibID":
        result=get_year(key)
        if result != -9999 :
          backlinks[key]+=backlink+";;"

          #DEBUG
          if key=="Amadon2008a":
            print(" Found key Amadon2008a, backlink:",backlink)
          #ENDDEBUG
          return '<a href="../../bibliography/generated_files/bibliography.html#%s">[%s]</a>' %(key,text)

    return '<a href="#">[[FAKE LINK:'+key+']]</a>'

  p=re.compile("\\[\\[([a-zA-Z0-9_ */<>.|:+#]*)\\]\\]")
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

def assemble_html(origin_yml_files,suppl_components,dir_name,root_filname,allowed_link_seeds,backlinks):
  """ Use the list of dictionaries "origin_yml_files" as well as the
      supplementary components "suppl_components", to produce html files,
      situated in dir_name (e.g. "tutorial") directories dir_name+"/generated_files".
      The root name of the files is "root_filname" (e.g. "lesson").
      The complementary yml information (intro, body, and possibly secX) for each file
      is situated in dir_name+"/origin_files".
      When there is a list of sections, a table of content is constituted automatically.
      The dictionary allowed_link_seeds allows one to set up the links to relevant keywords.
      The backlinks are accumulated, to be mentioned in the bibliography.html file.
      WARNING : not all files are assembled using this function ! In particular, the "topics" files are assembled in the main code ...
  """

  # Store the default informations
  for i, origin_yml in enumerate(origin_yml_files):
    if origin_yml.name.strip()=="default":
      origin_yml_default=origin_yml

  # Generate each "normal" file : assemble the content, apply global transformations, then write.
  list_names=[]
  dic_subtitles={}
  for i, origin_yml in enumerate(origin_yml_files):
    name = origin_yml.name
    if name=="default":
      continue
    list_names.append(name)
    dic_subtitles[name]=origin_yml.keyword+' - '+origin_yml.subtitle 
  
    if root_filname != "":
      full_filname=root_filname+"_"+name
    else:
      full_filname=name

    # Try to complete the information by reading a yml file
    path_ymlfile="%s/origin_files/%s.yml" %(dir_name,full_filname)
    doc_yml={}
    if os.path.isfile(path_ymlfile):
      print("Read "+path_ymlfile+" to build "+full_filname+".html ... ",end="")
      doc_yml=read_yaml(path_ymlfile)
 
    # If there are some sections in this yml file, constitute the body of the file using these sections,
    # and also make a table of content
    labels=[]
    for j in doc_yml.keys(): 
      if "sec" in j[:3]:
        labels.append(j[3:])
    secs_html=""
    if len(labels)!=0:
      labels.sort() 
      secs_html="\n <ul> \n"
      # Table of content
      for label in labels:
         secj="sec"+label
         secs_html+='  <li><a href="#%s">%s.</a> %s</li>\n' %(label,label,doc_yml[secj]["title"])
      secs_html+="\n </ul> \n <hr>"
      # Body
      for label in labels:
         secj="sec"+label
         sec_html='<p><a name="%s"> </a><br>' %(label)
         sec_html+='<h3><b>%s. %s</b></h3>\n <p>' %(label,doc_yml[secj]["title"])
         sec_html+=doc_yml[secj]["body"]
         full_filname=origin_yml.name
         if root_filname != "":
           full_filname=root_filname+"_"+origin_yml.name
         backlink=' &nbsp; <a href="../../%s/generated_files/%s.html#%s">%s#%s</a> &nbsp; ' %(dir_name,full_filname,label,full_filname,label)
         sec_html = make_links(sec_html,None,allowed_link_seeds,backlinks,backlink)
         secs_html+=sec_html
 
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
              "content","body","sections",
              "links","end"]:

      item=""
      # Try to get the item from different sources
      if j in doc_yml.keys() and (not "sec" in j or j=="sections"):
        item+=doc_yml[j]
      elif j in suppl.keys() :
        item+=suppl[j]
      elif j =="sections":
        item+=secs_html
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

    rc=finalize_html(doc_html,origin_yml,dir_name,root_filname,allowed_link_seeds,backlinks) 

  # Generate a html list of files. WARNING : the "topics" list of files is NOT generated here ...
  list_names.sort()
  toc_all = '<a name="list"></a>'
  toc_all += '<h3><b> Alphabetical list of all files.</b></h3>'
  cur_let_all = 'A'
  toc_all += "<p>"+cur_let_all+".\n"

  for (ii,name) in enumerate(list_names):
    while not (name.startswith(cur_let_all.lower()) or name.startswith(cur_let_all.upper())):
      cur_let_all = chr(ord(cur_let_all)+1)
      toc_all = toc_all + "<p>"+cur_let_all+".\n"

    if root_filname != "":
      full_filname=root_filname+"_"+name
    #This is exceptional, should be soon removed
    elif name=="external" or name=="allvars":
      full_filname="varset_"+name
    else:
      full_filname=name

    toc_all = toc_all + '<br><a href="%s.html"/>%s</a> [%s] \n' %(full_filname,name,dic_subtitles[name])

  all_files_html=""
  spec={'users':'help files','tutorial':'lessons of the tutorial',
        'theory':'theory documents','input_variables':'varsets','bibliography':'generated files in the bibliography directory'}
  for j in ["header","title","subtitle","copyright","links","toc_all","links","end"]:
    if j == "toc_all":
      all_files_html += toc_all
    elif j == "subtitle":
      all_files_html += 'This document lists the %s of the ABINIT package.' %(spec[dir_name])
    else:
      all_files_html += getattr(origin_yml_default,j)

  root_filname_all="all_files"
  html_section=""
  rc=finalize_html(all_files_html,origin_yml_default,dir_name,root_filname_all,allowed_link_seeds,backlinks)

  return "Exit assemble_html"

################################################################################

def finalize_html(doc_html,origin_yml,dir_name,root_filname,allowed_link_seeds,backlinks):
  """ Final steps of the preparation of a html file : global operations and write
      The draft text of the file to be written is contained in doc_html
      Some substitutions must be done using information in origin_yml (.name, .howto, .keyword, .authors)
      The links are build using allowed_link_seeds.
      The backlinks are given back.
  """

  if root_filname == "all_files":
    full_filname=root_filname
  elif root_filname != "":
    full_filname=root_filname+"_"+origin_yml.name
  else:
    full_filname=origin_yml.name
  backlink=' &nbsp; <a href="../../%s/generated_files/%s.html">%s</a> &nbsp; ' %(dir_name,full_filname,full_filname)

  doc_html=doc_html.replace("__JS_PATH__","../../js_files/")

  if origin_yml.howto != "":
    doc_html=doc_html.replace("__HOWTO__",origin_yml.howto)

  if origin_yml.keyword != "":
    doc_html=doc_html.replace("__KEYWORD__",origin_yml.keyword)

  if origin_yml.authors != "":
    doc_html=doc_html.replace("__AUTHORS__",origin_yml.authors)

  doc_html = make_links(doc_html,None,allowed_link_seeds,backlinks,backlink)

  #This must be temporary !
  if full_filname=='allvars':
    full_filname='varset_allvars'
  if full_filname=='external':
    full_filname='varset_external'

  #Write the finalized html file.
  path_html = "%s/generated_files/%s.html" %(dir_name,full_filname)
  file_html = open(path_html,'w')
  file_html.write(doc_html)
  file_html.close()
  print("File %s written ..."%path_html)

  return "Exit finalize_html"

################################################################################
