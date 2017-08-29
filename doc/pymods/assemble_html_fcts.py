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

def format_backlinks(backlinks_str):
  """ Post-process the list of backlinks accumulated for a particular reference 
  """
  backlinks_list=backlinks_str.split(";;")
  backlinks_newlist=[]
  if len(backlinks_list)!=0:
    list_stripped=[]
    for link in backlinks_list:
      stripped=link.strip()
      if stripped!="":
        list_stripped.append(stripped)
    if len(list_stripped)!=0:
      set_stripped=set(list_stripped)
      if len(set_stripped)!=0:
        backlinks_newlist=list(set_stripped)
        backlinks_newlist.sort()
      else:
        backlinks_newlist=[]
    else:
      backlinks_newlist=[]
  backlinks_formatted=""
  if len(backlinks_newlist)!=0:
    for ilink in range(len(backlinks_newlist)-1):
      backlinks_formatted+=backlinks_newlist[ilink]+", "
    backlinks_formatted+=backlinks_newlist[len(backlinks_newlist)-1]+". "
  return backlinks_formatted

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
    #Explicitly : namespace:key#section|dokutext
    #First, the namespace
    if ':' in dokukey:
      p1234=dokukey.split(':',1)
      namespace=p1234[0].strip()
      p234=p1234[1].strip()
    else:
      namespace=""
      p234=dokukey
    #Then, the dokutext
    if '|' in p234:
      p234_split=p234.split('|',1) 
      p23=p234_split[0].strip()
      dokutext=p234_split[1].strip()
    else:
      p23=p234
      dokutext=""
    #Finally, the key (often a filename, but not for input variables) and section
    if '#' in p23:
      p23_split=p23.split('#',1)
      key=p23_split[0].strip()
      section=p23_split[1].strip()
    else:
      key=p23
      section=""

    #Detect cases of external namespace
    external_namespace=0
    if namespace!="" and namespace in ["http","https","ftp","file"]:
      external_namespace=1

    #Prepare the webtext in the simple cases, note the namespace is echoed only when it is external
    webtext=dokutext
    if webtext=="":
      if external_namespace==1:
        webtext+=namespace+":"
      if key!="":
        webtext+=p23
      else:
        webtext+=section

    #Finalize the cases of external links
    if external_namespace==1:
      return '<a href="%s:%s">%s</a>' %(namespace,p23,webtext)
    if namespace=="" and key[:4]=="www.":
      return '<a href="http://%s">%s</a>' %(p23,webtext)

    #Treat the internal links
    if namespace=="":
      linkseed=key
      if key=="":
        #This is the own file, no need to establish backlinks or further processing, the section is the reference
        return '<a href="#%s">%s</a>' %(section,webtext)
    elif namespace in ["aim","anaddb","optic"]:
      linkseed=key+"@"+namespace
    else:
      linkseed=namespace+'_'+key

    #The allowed namespaces are:
    dic_namespaces={"aim":"input_variables/generated_files",
                    "anaddb":"input_variables/generated_files",
                    "bib":"biblio/generated_files",
                    "help":"users/generated_files",
                    "lesson":"tutorial/generated_files",
                    "optic":"input_variables/generated_files",
                    "theorydoc":"theory/generated_files",
                    "topic":"topics/generated_files",
                    "varset":"input_variables/generated_files"}

    #Actually for the internal links, make the selection on the linkseed at present ... this should be changed ...
    #Might be changed, later ...
    if linkseed in allowed_link_seeds.keys():
      value=allowed_link_seeds[linkseed]

      #DEBUG
      #if "help:respfn" in dokukey:
      #  print(" ")
      #  print(" dokukey:",dokukey)
      #  print(" value:",value)
      #  print(" ")
      #ENDDEBUG

      #Treat first the allowed namespaces
      if value in dic_namespaces.keys():
        dir=dic_namespaces[value]

        #Specific formatting treatment
        if value=="help" and namespace=="":
          webtext=key[5:]+' help file'
        
        #Set up backlink for topic (might be generalized to all ?)
        if value=="topic":
          backlinks[linkseed]+=backlink+";;"

        return '<a href="../../%s/%s.html#%s">%s</a>' %(dir,linkseed,section,webtext)

      #Treat everything else
      elif "input_variable in " in value:
        # This is a link to an input variable
        filename=value[18:]
        return '<a href="../../input_variables/generated_files/varset_%s.html#%s">%s</a>' %(filename,linkseed,webtext)
      elif value=="characteristic":
        return '<a href="../../users/generated_files/help_abinit.html#%s">%s</a>' %(key,webtext)
      elif value=="in_tests":
        return '<a href="../../%s">&#126;abinit/%s</a>' %(key,key)
      elif value=="bibID":
        result=get_year(key)
        if result != -9999 :
          backlinks[key]+=backlink+";;"
          return '<a href="../../biblio/generated_files/bib_biblio.html#%s">[%s]</a>' %(key,webtext)

    return '<a href="#">[[FAKE LINK:'+dokukey+']]</a>'

  p=re.compile("\\[\\[([a-zA-Z0-9_ ,\(\)\"\'?=*/<>.|:+#@-]*)\\]\\]")
  if text is None:
    return ""
  new_text=p.sub(replace_link,text)

  return new_text

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
      The backlinks are accumulated, to be mentioned in the bib_biblio.html file.
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
    # and also make a table of content. Allows (only) two levels for the table of content.
    labels=[]
    for j in doc_yml.keys(): 
      if "sec" in j[:3]:
        labels.append(j[3:])
    secs_html=""
    if len(labels)!=0:
      #Trick (not perfect ...) to sort numbers and digits together. Will not work with strings longer than 9 digits.
      labels=sorted(labels, key= lambda item: str(len(item))+item)
      secs_html="\n <ul> \n"
      # Table of content
      for label in labels:
         secj="sec"+label
         secs_html+='  <li>%s. <a href="#%s">%s</a></li>\n' %(label,doc_yml[secj]["tag"],doc_yml[secj]["title"])
         #Treat one level of subsections
         sublabels=[]
         for subj in doc_yml[secj].keys():
           if "sec" in subj[:3]:
             sublabels.append(subj[3:])
         if len(sublabels)!=0:
           sublabels=sorted(sublabels, key= lambda item: str(len(item))+item)
           secs_html+="\n   <ul> \n"
           for sublabel in sublabels:
             subsecj="sec"+sublabel
             secs_html+='    <li>%s. <a href="#%s">%s</a></li>\n' %(sublabel,doc_yml[secj][subsecj]["tag"],doc_yml[secj][subsecj]["title"])
           secs_html+="\n   </ul> \n"
      secs_html+="\n </ul> \n <hr>"
      # Body
      for label in labels:
         secj="sec"+label
         sec_html='<br><a name="%s"> </a>\n' %(doc_yml[secj]["tag"])
         sec_html+='<a name="%s"> </a>\n' %(label)
         sec_html+='<h3><b>%s. %s</b></h3>\n <p>' %(label,doc_yml[secj]["title"])
         if "body" in doc_yml[secj].keys():
           sec_html+=doc_yml[secj]["body"]
           full_filname=origin_yml.name
           if root_filname != "":
             full_filname=root_filname+"_"+origin_yml.name
           backlink=' &nbsp; <a href="../../%s/generated_files/%s.html#%s">%s#%s</a> ' %(dir_name,full_filname,label,full_filname,label)
           sec_html = make_links(sec_html,None,allowed_link_seeds,backlinks,backlink)
         #Treat one level of subsections
         sublabels=[]
         for subj in doc_yml[secj].keys():
           if "sec" in subj[:3]:
             sublabels.append(subj[3:])
         if len(sublabels)!=0:
           sublabels=sorted(sublabels, key= lambda item: str(len(item))+item)
           for sublabel in sublabels:
             subsecj="sec"+sublabel
             sec_html+='<br><a name="%s"> </a>\n' %(doc_yml[secj][subsecj]["tag"])
             sec_html+='<a name="%s"> </a>\n' %(sublabel)
             sec_html+='<h4><b>%s. %s</b></h4>\n <p>' %(sublabel,doc_yml[secj][subsecj]["title"])
             sec_html+=doc_yml[secj][subsecj]["body"]
             full_filname=origin_yml.name
             if root_filname != "":
               full_filname=root_filname+"_"+origin_yml.name
             backlink=' &nbsp; <a href="../../%s/generated_files/%s.html#%s">%s#%s</a> ' %(dir_name,full_filname,sublabel,full_filname,sublabel)
             sec_html = make_links(sec_html,None,allowed_link_seeds,backlinks,backlink)
         sec_html+="<br><br><a href=#top>Go to the top</a>\n<hr>\n"
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
        'theory':'theory documents','input_variables':'varsets','biblio':'generated files in the biblio directory'}
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
  backlink=' &nbsp; <a href="../../%s/generated_files/%s.html">%s</a> ' %(dir_name,full_filname,full_filname)

  doc_html=doc_html.replace("__JS_PATH__","../../js_files/")

  if origin_yml.howto != "":
    doc_html=doc_html.replace("__HOWTO__",origin_yml.howto)

  if origin_yml.keyword != "":
    doc_html=doc_html.replace("__KEYWORD__",origin_yml.keyword)

  if origin_yml.authors != "":
    doc_html=doc_html.replace("__AUTHORS__",origin_yml.authors)

  doc_html = make_links(doc_html,None,allowed_link_seeds,backlinks,backlink)

  #Write the finalized html file.
  path_html = "%s/generated_files/%s.html" %(dir_name,full_filname)
  file_html = open(path_html,'w')
  file_html.write("<!--                                                                                     -->\n")
  file_html.write("<!-- This file has been produced by generate_doc.py using different .yml  or .bib files. -->\n")
  file_html.write("<!-- It is useless to modify it. Modify the related .yml  or .bib files instead.         -->\n")
  file_html.write("<!--                                                                                     -->\n")
  try:
    file_html.write(doc_html)
  except:
    #print("\n\n... ERROR ...\n[Complement of information from generate_doc.py]")
    #print(" doc_html[]:",doc_html[3800:3850],"\n")
    raise

  file_html.close()
  print("File %s written ..."%path_html)

  return "Exit finalize_html"

################################################################################
