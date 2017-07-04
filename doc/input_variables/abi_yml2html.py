#! /usr/bin/env python 
from __future__ import print_function

import sys
import os
import yaml
import re
import argparse
from variables import *

debug = 0
make_topics_visible=1

# Path relative from HTML files
js_path = "../"
users_path = "../../users/"

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

def make_links(text,cur_abivarname,variables,characteristics,specials):

  def replace_link(mymatch):
    abivarname = mymatch.group()[2:-2]
    if abivarname == cur_abivarname:
      return "<b>"+cur_abivarname+"</b>"
    elif abivarname in variables.keys():
      section = variables[abivarname]
      return "<a href=\""+section+".html#"+abivarname+"\">"+abivarname+"</a>"
    elif abivarname in characteristics:
      return "<a href=\""+users_path+"abinit_help.html#"+str.replace(abivarname.lower()," ","_")+"\">"+abivarname+"</a>"
    elif abivarname in specials:
      return "<a href=\"specials.html#"+abivarname+"\">"+abivarname+"</a>"
    else:
      return "<a href=\"#\">[[FAKE LINK:"+abivarname+"]]</a>"
    return mymatch.group()

  p=re.compile("\\[\\[([a-zA-Z0-9_ */<>]*)\\]\\]")
  if text is None:
    return ""
  new_text=p.sub(replace_link,text)

  return new_text

################################################################################

def read_yaml_file(ymlfile):
  """ Read the file 'ymlfile', containing yml data, and store all such data in the returned object"""

  if ymlfile== "abinit_vars.yml":
    print("Will use abinit_vars.yml as database input file for the input variables and their characteristics ...")
  elif ymlfile== "yml_files/list_of_topics.yml":
    print("Will use yml_files/list_of_topics.yml as database input file for the list of topics ...")
  elif ymlfile== "topics.yml":
    print("Will use topics.yml as database input file for the list of topics ...")
  elif ymlfile== "sections.yml":
    print("Will use sections.yml as database input file for the list of sections ...")
  elif ymlfile== "yml_files/tests_dirs.yml":
    print("Will use yml_files/tests_dirs.yml as database input file for the list of directories in which automatic test input files are present ...")
  elif ymlfile== "html_automatically_generated/topics_in_tests.yml":
    print("Generated file named html_automatically_generated/topics_in_tests.yml, to contain the list of automatic test input files relevant for each topic ...")

  with open(ymlfile, 'r') as f:
    ymlstructure = yaml.load(f);

  return ymlstructure

################################################################################
################################################################################

# Parsing section

################################################################################
 
variables=read_yaml_file("abinit_vars.yml")
list_of_topics=read_yaml_file("yml_files/list_of_topics.yml")
topics=read_yaml_file("topics.yml")
sections=read_yaml_file("sections.yml")
tests_dirs=read_yaml_file("yml_files/tests_dirs.yml")

################################################################################
# Parse the ABINIT input files, in order to find the possible topics to which they are linked -> topics_in_tests

try :
  rm_cmd = "rm html_automatically_generated/topics_in_tests.yml"
  retcode = os.system(rm_cmd)
except :
  if debug==1 :
    print("rm html_automatically_generated/topics_in_tests.yml failed")
    print("the file was likely non existent")

for tests_dir in tests_dirs :
  grep_cmd = "grep topics tests/%s/Input/*.in > html_automatically_generated/topics_in_tests.txt" % (tests_dir)
  retcode = os.system(grep_cmd)
  if retcode == 0 :
    sed_cmd = "sed -e 's/^/- /' html_automatically_generated/topics_in_tests.txt >> html_automatically_generated/topics_in_tests.yml "
    retcode = os.system(sed_cmd)

topics_in_tests=read_yaml_file("html_automatically_generated/topics_in_tests.yml")
if debug==1 :
  print(" topics_in_tests :")
  print(topics_in_tests)

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
  cur_content = "<br><font id=\"title\"><a name=\""+speckey+"\">"+speckey+"</a></font>\n"
  cur_content += "<br><font id=\"text\">\n"
  cur_content += "<p>\n"+doku2html(make_links(specval,speckey,list_all_vars,list_chars,cur_specials))+"\n"
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
      cur_content += "<br><font id=\"characteristic\">Characteristic: "+make_links(chars,var.abivarname,list_all_vars,list_chars,cur_specials)+"</font>\n"
    else:
      cur_content += "<br><font id=\"characteristic\">Characteristic: </font>\n"
    # Topics
    try:
      if var.topics is not None and make_topics_visible==1 :
        cur_content += "<br><font id=\"characteristic\">Mentioned in \"How to\": "
        vartopics=var.topics
        topics_name_class = vartopics.split(',')
        for i, topic_name_class in enumerate(topics_name_class):
          name_class = topic_name_class.split('_')
          cur_content += "<a href=\""+name_class[0]+".html\">"+name_class[0]+"</a> "
        cur_content += "</font>\n"
    except:
      if debug==1 :
        print(" No topic_class for abivarname "+var.abivarname)
    # Variable type, including dimensions
    cur_content += "<br><font id=\"vartype\">Variable type: "+var.vartype
    if var.dimensions is not None:
      cur_content += make_links(format_dimensions(var.dimensions),var.abivarname,list_all_vars,list_chars,cur_specials)
    if var.commentdims is not None and var.commentdims != "":
      cur_content += " (Comment: "+make_links(var.commentdims,var.abivarname,list_all_vars,list_chars,cur_specials)+")"
    cur_content += "</font>\n" 
    # Default
    cur_content += "<br><font id=\"default\">"+make_links(format_default(var.defaultval),var.abivarname,list_all_vars,list_chars,cur_specials)
    if var.commentdefault is not None and var.commentdefault != "":
      cur_content += " (Comment: "+make_links(var.commentdefault,var.abivarname,list_all_vars,list_chars,cur_specials)+")"
    cur_content += "</font>\n" 
    # Requires
    if var.requires is not None and var.requires != "":
      cur_content += "<br><br><font id=\"requires\">\nOnly relevant if "+doku2html(make_links(var.requires,var.abivarname,list_all_vars,list_chars,cur_specials))+"\n</font>\n"
    # Excludes
    if var.excludes is not None and var.excludes != "":
      cur_content += "<br><br><font id=\"excludes\">\nThe use of this variable forbids the use of "+doku2html(make_links(var.excludes,var.abivarname,list_all_vars,list_chars,cur_specials))+"\n</font>\n"
    # Text
    cur_content += "<br><font id=\"text\">\n"
    cur_content += "<p>\n"+doku2html(make_links(var.text,var.abivarname,list_all_vars,list_chars,cur_specials))+"\n"
    # End the section for one variable
    cur_content += "</font>\n\n"
    cur_content += "<br><br><br><br><a href=#top>Go to the top</a>\n"
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
  sectionhtml = doku2html(make_links(sectionhtml,None,list_all_vars,list_chars,cur_specials))

  #Write the finalized html file.
  file_cur = 'html_automatically_generated/'+section+'.html'
  f_cur = open(file_cur,'w')
  f_cur.write(sectionhtml)
  f_cur.write("\n")
  f_cur.close()
  print("File "+section+".html has been written ...")

################################################################################
################################################################################

# Generate the different "topic" files

################################################################################
# Constitute the section "Related input variables" for all topic files. 
# This section in input variables is stored, for each topic_name, in topic_invars[topic_name]

topic_invars = dict()
found = dict()

for i, topic in enumerate(topics):
  topic_invars[topic.topic_name] = ""

for (tclasskey, tclassval) in list_topics_class:

  if debug == 1:
    print("\nWork on "+tclasskey+"\n")
  for i, topic in enumerate(topics):
    found[topic.topic_name] = 0

  for i, var in enumerate(variables):
    try:
      if var.topics is not None:
        topics_name_class = var.topics.split(',')
        for i, topic_name_class in enumerate(topics_name_class):
          name_class = topic_name_class.split('_')
          if tclasskey==name_class[1].strip() :
            topic_name=name_class[0].strip()
            if found[topic_name]==0 :
              topic_invars[topic_name] += "<p>"+tclassval+"<p>"
              found[topic_name] = 1
            abivarname=var.abivarname
            if var.characteristics is not None and '[[INTERNAL_ONLY]]' in var.characteristics:
              abivarname = '%'+abivarname
            topic_invars[topic_name] += "... <a href=\""+var.section+".html#"+var.abivarname+"\">"+abivarname+"</a>   "
            topic_invars[topic_name] += "["+var.mnemonics+"]<br>\n"
    except:
      if debug==1 :
       print(" No topics for abivarname "+var.abivarname) 

################################################################################
# Constitute the section 4 "Selected input files" for all topic files.
# This section on input files in topics  is stored, for each topic_name, in topic_infiles[topic_name]

topic_infiles = dict()
for i, var in enumerate(topics):
  topic_infiles[var.topic_name] = ""
 
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
  tests=inputs_for_topic[topic_name]
# Constitute a dictionary for each directory of tests
  dir = dict()
  for test in tests:
    file_split=test.split('/')
    dirname=file_split[1] 
    testname=file_split[3] 
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
# Assemble the "topic" files 

print("NEW : Will use file yml_files/default_topic.yml as default for all topic files ... ")
default_topic_yml=read_yaml_file("yml_files/default_topic.yml")
default_topic=default_topic_yml[0]

# For each "topic" file
for topic_name in list_of_topics:
  f_newtopic="yml_files/topic_"+topic_name+".yml"
  print("NEW : Will use file "+f_newtopic+" to initiate the topic "+topic_name+" ... ")
  newtopic_yml=read_yaml_file("yml_files/topic_"+topic_name+".yml")
  newtopic=newtopic_yml[0]

  #Generate the table of content
  item_toc=0
  item_list=[]
  title={ "introduction":"Introduction." , "examples":"Example(s)", "tutorials":"Related tutorials." , "input_variables":"Related input variables." , "input_files":"Selected input files."}
  sec_number={}
  toc=" <h3><b>Table of content: </b></h3> \n <ul> "
  for j in ["introduction","examples","tutorials","input_variables","input_files"] :
    sec_number[j]="0"
    try :
      extract_j=getattr(newtopic,j).strip()
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
      extract_j=getattr(newtopic,j).strip()
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
  topic_html=topic_html.replace("__HOWTO__",newtopic.howto)
  topic_html=topic_html.replace("__KEYWORD__",newtopic.keyword)
  topic_html = doku2html(make_links(topic_html,None,list_all_vars,list_chars,cur_specials))

  # Open, write and close the file
  file_topic = 'html_automatically_generated/topic_'+topic_name+'.html'
  f_newtopic = open(file_topic,'w')
  f_newtopic.write(topic_html)
  f_newtopic.close()
  print("File topic_"+topic_name+".html has been written ...")

################################################################################
# Generate the file with the list of names of different "topic" files

#Constitute first the table of content
toc = '<a name="list"></a>'
toc += '<h3><b> Alphabetical list of all "How to ?" documentation topics.</b></h3>'

cur_let = 'A'
toc += "<p>"+cur_let+".&nbsp;\n"
for i, topic in enumerate(topics):
  while not (topic.topic_name.startswith(cur_let.lower()) or topic.topic_name.startswith(cur_let.upper())):
    cur_let = chr(ord(cur_let)+1)
    toc = toc + "<p>"+cur_let+".\n"
  topic_name = topic.topic_name
  toc = toc + "<br><a href=\"topic_"+ topic_name + ".html\">" + topic_name + "</a> [How to "+topic.howto+"] &nbsp;&nbsp;\n"

#Generate a first version of the html file, in the order "header" ... up to the "end"
#Take the info from the section "default" if there is no information on the specific section provided in the yml file.
all_topics_html=""
for j in ["header","title","subtitle","copyright","links","toc","links","end"]:
  if j == "toc":
    all_topics_html += toc
  elif j == "subtitle":
    all_topics_html += "<h2>Complete list.</h2><hr>"
    all_topics_html += 'This document lists the names of all "How to ?" documentation topics for the abinit package.'
    all_topics_html += '<script type="text/javascript" src="../generic_advice.js"> </script>'
  else:
    all_topics_html += getattr(default_topic,j)

#Global operations on the tentative html file.
all_topics_html=all_topics_html.replace("__JS_PATH__",js_path)
all_topics_html=all_topics_html.replace("__KEYWORD__",default_topic.keyword)
all_topics_html = doku2html(make_links(all_topics_html,None,list_all_vars,list_chars,cur_specials))

# Open, write and close the file
file_html = 'html_automatically_generated/all_topics.html'
f_html = open(file_html,'w')
f_html.write(all_topics_html)
f_html.close()
print("File all_topics.html has been written ...")
print("Work done !")

################################################################################
# Assemble the "topic" files : secs 1, 2, 3 and then possibly 4 and 5

# For each "topic" file 
for topic_name, content in topic_invars.items():
  file_template = 'html_template/temp_'+topic_name+'.html'
  print("Will use file named temp_"+topic_name+", as template for "+topic_name+".html... ", end='')
  try:
    with open(file_template) as f:
      header_varX = f.read()
  except:
    print("Tried to open the file temp_"+topic_name+" , but failed.")

  file_topic = 'html_automatically_generated/'+topic_name+'.html'
  topic_header_varX = header_varX.replace("__JS_PATH__",js_path)
  if topic_name not in "GS_introduction":
    topic_header_varX += "<script type=\"text/javascript\" src=\""+js_path+"related_input_variables.js\"> </script>\n\n"

  f_topic = open(file_topic,'w')

  topic_header_varX2 = doku2html(make_links(topic_header_varX,None,list_all_vars,list_chars,cur_specials))

  f_topic.write(topic_header_varX2)

# Write Sec. 3 
#  f_topic.write("\n\n<p>&nbsp; \n<HR ALIGN=left> \n<p> <a name=\"3\">&nbsp;</a>\n<h3><b> 3. Related input variables.</b></h3>\n\n\n")
  f_topic.write(topic_invars[topic_name])

# Write Sec. 4
  f_topic.write("\n\n<p>&nbsp; \n<HR ALIGN=left> \n<p> <a name=\"4\">&nbsp;</a>\n<h3><b> 4. Selected input files.</b></h3>\n\n\n")
  f_topic.write("The user can find some related example input files in the ABINIT package in the directory /tests, or on the Web:\n")
  f_topic.write(topic_infiles[topic_name])

# Write final lines
  content = "<br>"
  content += "<script type=\"text/javascript\" src=\""+js_path+"list_internal_links_incl_topics.js\"> </script>\n\n"
  content += "</body>\n"
  content += "</html>"
  content += "\n"
  f_topic.write(content)
  f_topic.close()
  print("File "+topic_name+".html has been written ...")

################################################################################
################################################################################

# Generate the file with the list of names of different "topic" files

#Generate the files that document all the variables (all such files : var* as well as all and special).
# Parse the header of alltopics file and also replace the JS_PATH.

with open('html_template/temp_alltopics.html') as f:
    header_alltopics = f.read()
print("Will use file named temp_alltopics.html, as template for alltopics.html... ", end='')

header_alltopics = header_alltopics.replace("__JS_PATH__",js_path)

# Initialize the alltopics file output
toutput = ''
toutput = toutput + header_alltopics + "<br />\n"
cur_let = 'A'
toutput = toutput + "<p>"+cur_let+".&nbsp;\n"

# DEBUG Test the content of topics
# for topic in topics:
#  print("topic: "+topic.topic_name)

if debug==1 :
  print("Will enter loop on topics")

for i, topic in enumerate(topics):
  if debug==1 :
    print("topic: "+topic.topic_name)
    print("cur_let: "+cur_let)
  while not (topic.topic_name.startswith(cur_let.lower()) or topic.topic_name.startswith(cur_let.upper())):
    cur_let = chr(ord(cur_let)+1)
    toutput = toutput + "<p>"+cur_let+".\n"
  topic_name = topic.topic_name
  toutput = toutput + "<br><a href=\""+ topic_name + ".html\">" + topic_name + "</a> [How to "+topic.howto+"] &nbsp;&nbsp;\n"

################################################################################
# Alltopics file : complete the content, then write the file and close it
toutput += "<script type=\"text/javascript\" src=\""+js_path+"list_internal_links_incl_topics.js\"> </script>\n\n"
toutput += "</body>\n"
toutput += "</html>"

file_html = 'html_automatically_generated/alltopics.html'
f_html = open(file_html,'w')
f_html.write(toutput)
f_html.close()
print("File alltopics.html has been written ...")
print("Work done !")

################################################################################
