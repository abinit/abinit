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
# Definitions

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

def doku2html(text):

  def replace_link(mymatch):
    abivarname = mymatch.group()[2:-2]
    return "<b>"+abivarname+"</b>"

  p = re.compile("\*\*([a-zA-Z0-9_ */<>]*)\*\*")
  text2 = p.sub(replace_link,text)

  return text2

def format_default(defaultval):

  if defaultval is None:
    s = 'No default'
  else:
    s = "Default is "+str(defaultval)
 
  return s

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
  #try:
  new_text=p.sub(replace_link,text)
  #except:
  #  print("Trying to compile :"+str(text))


  return new_text

################################################################################
################################################################################

# Parsing section

################################################################################
# Parse the abinit_vars.yml file -> variables

file='abinit_vars.yml'
print("Will use "+str(file)+" as database input file for the input variables and their characteristics ...")

parser = argparse.ArgumentParser(description='Tool for eigenvalue analysis')
parser.add_argument('-f','--file',help='YML file to be read')
args = parser.parse_args()
args_dict = vars(args)
if args_dict['file']:
  file = args_dict['file']

with open(file, 'r') as f:
    variables = yaml.load(f);


################################################################################
# Parse the topics.yml file -> topics

file='topics.yml'
print("Will use "+str(file)+" as database input file for the list of topics ...")

parser = argparse.ArgumentParser(description='Tool for eigenvalue analysis')
parser.add_argument('-f','--file',help='YML file to be read')
args = parser.parse_args()
args_dict = vars(args)
if args_dict['file']:
  file = args_dict['file']

with open(file, 'r') as f:
    topics = yaml.load(f);

################################################################################
# Parse the tests_dirs.yml file -> tests_dirs

file='tests_dirs.yml'
print("Will use "+str(file)+" as database input file for the list of directories in which automatic test input files are present ...")

parser = argparse.ArgumentParser(description='Tool for eigenvalue analysis')
parser.add_argument('-f','--file',help='YML file to be read')
args = parser.parse_args()
args_dict = vars(args)
if args_dict['file']:
  file = args_dict['file']

with open(file, 'r') as f:
    tests_dirs = yaml.load(f);

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

file='html_automatically_generated/topics_in_tests.yml'
print("Generated file named "+str(file)+", to contain the list of automatic test input files relevant for each topic ...")

parser = argparse.ArgumentParser(description='Tool for eigenvalue analysis')
parser.add_argument('-f','--file',help='YML file to be read')
args = parser.parse_args()
args_dict = vars(args)
if args_dict['file']:
  file = args_dict['file']

with open(file, 'r') as f:
  topics_in_tests = yaml.load(f);

if debug==1 :
  print(" topics_in_tests :")
  print(topics_in_tests)

################################################################################
# Parse the headers of allvariables and special variables files also replace the JS_PATH.

with open('html_template/temp_allvariables.html') as f:
    header_all = f.read()
print("Will use file named temp_allvariables.html as template for allvariables.html...")

with open('html_template/temp_specials.html') as f:
    header_specials = f.read()
print("Will use file named temp_specials.html as template for specials.html...")

header_all = header_all.replace("__JS_PATH__",js_path)

all_contents = dict()
all_vars = dict()

# Create a dictionary that gives the section for each abivarname

list_all_vars = dict()

for var in variables:
  list_all_vars[var.abivarname] = var.section

#
cur_specials = []
for (specialkey,specialval) in list_specials:
  cur_specials.append(specialkey)

################################################################################
# Constitute the body of information for all variables, stored for the appropriate section in all_contents[section]

# Initialize the all variable file output
output = ''
output = output + header_all + "<br />\n"
cur_let = 'A'
output = output + "<p>"+cur_let+".&nbsp;\n"

for i, var in enumerate(variables):
  if debug==1 :
    print(var)
  while not var.abivarname.startswith(cur_let.lower()):
    cur_let = chr(ord(cur_let)+1)
    output = output + "<p>"+cur_let+".\n"
  abivarname = var.abivarname
  if var.characteristics is not None and '[[INTERNAL_ONLY]]' in var.characteristics:
    abivarname = '%'+abivarname
  output = output + "<a href=\""+var.section+".html#"+var.abivarname+"\">"+abivarname+"</a>&nbsp;&nbsp;\n" 
  section = var.section
  if section not in all_contents.keys():
    all_contents[section] = "<br><br><br><br><hr>\n"
    all_vars[section] = []

  all_vars[section].append([var.abivarname,var.mnemonics])

  # Constitute the body of information related to one input variable
  cur_content = ""

  try:
    cur_content += "<br><font id=\"title\"><a name=\""+var.abivarname+"\">"+var.abivarname+"</a></font>\n"
    cur_content += "<br><font id=\"mnemonics\">Mnemonics: "+var.mnemonics+"</font>\n"
    if var.characteristics is not None:
      chars = ""
      for chs in var.characteristics:
        chars += chs+", "
      chars = chars[:-2]
      cur_content += "<br><font id=\"characteristic\">Characteristic: "+make_links(chars,var.abivarname,list_all_vars,list_chars,cur_specials)+"</font>\n"
    else:
      cur_content += "<br><font id=\"characteristic\">Characteristic: </font>\n"
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
    cur_content += "<br><font id=\"vartype\">Variable type: "+var.vartype
    if var.dimensions is not None:
      cur_content += make_links(format_dimensions(var.dimensions),var.abivarname,list_all_vars,list_chars,cur_specials)
    if var.commentdims is not None and var.commentdims != "":
      cur_content += " (Comment: "+make_links(var.commentdims,var.abivarname,list_all_vars,list_chars,cur_specials)+")"
    cur_content += "</font>\n" 
    cur_content += "<br><font id=\"default\">"+make_links(format_default(var.defaultval),var.abivarname,list_all_vars,list_chars,cur_specials)
    if var.commentdefault is not None and var.commentdefault != "":
      cur_content += " (Comment: "+make_links(var.commentdefault,var.abivarname,list_all_vars,list_chars,cur_specials)+")"
    cur_content += "</font>\n" 
    if var.requires is not None and var.requires != "":
      cur_content += "<br><br><font id=\"requires\">\nOnly relevant if "+doku2html(make_links(var.requires,var.abivarname,list_all_vars,list_chars,cur_specials))+"\n</font>\n"
    if var.excludes is not None and var.excludes != "":
      cur_content += "<br><br><font id=\"excludes\">\nThe use of this variable forbids the use of "+doku2html(make_links(var.excludes,var.abivarname,list_all_vars,list_chars,cur_specials))+"\n</font>\n"
    cur_content += "<br><font id=\"text\">\n"
    cur_content += "<p>\n"+doku2html(make_links(var.text,var.abivarname,list_all_vars,list_chars,cur_specials))+"\n"
    cur_content += "</font>\n\n"
    cur_content += "<br><br><br><br><a href=#top>Go to the top</a>\n"
    cur_content += "<B> | </B><a href=\"allvariables.html#top\">Complete list of input variables</a><hr>\n"
    all_contents[section] = all_contents[section] + cur_content + "\n\n"
  except AttributeError as e:
    print(e)
    print('For variable : ',abivarname)

################################################################################
# Generate the files that document all the variables.

# For each "normal" section file : generate the header, generate the alphabetical list, write these,
# then complete the content (that was previously gathered), then write the file and close it
for section, content in all_contents.items():
  file_cur = 'html_automatically_generated/'+section+'.html'
  f_cur = open(file_cur,'w')

  with open('html_template/temp_'+section+'.html') as f:
    header_varX = f.read()
  print("Will use file named temp_"+section+", as template for "+section+".html... ", end='')

  cur_header_varX = header_varX.replace("__JS_PATH__",js_path)
 
  f_cur.write(cur_header_varX)
  cur_let = 'A'
  f_cur.write(" <br>"+cur_let+".\n")
  for abivarname,defi in all_vars[section]:
    #curlink = "<br /><a href=\"#"+abivarname+"\">"+abivarname+"</a> : "+defi+"\n"
    while not abivarname.startswith(cur_let.lower()):
      cur_let = chr(ord(cur_let)+1)
      f_cur.write(" <br>"+cur_let+".\n")
    curlink = " <a href=\"#"+abivarname+"\">"+abivarname+"</a>&nbsp;&nbsp;\n"
    f_cur.write(curlink)
    #f_cur.write("<br /><br />\n")
  f_cur.write("\n")
  if make_topics_visible==0:
    content += "<script type=\"text/javascript\" src=\""+js_path+"list_internal_links.js\"> </script>\n\n"
  else:
    content += "<script type=\"text/javascript\" src=\""+js_path+"list_internal_links_incl_topics.js\"> </script>\n\n"
  content += "</body>\n"
  content += "</html>"
  f_cur.write(content)
  f_cur.write("\n")
  f_cur.close()
  print("File "+section+".html has been written ...")

# Special file : complete the content, then write the file and close it
file_specials = 'html_automatically_generated/specials.html'
f_sp = open(file_specials,'w')
header_specials = header_specials.replace("__JS_PATH__",js_path)
f_sp.write(header_specials)
cur_let = 'A'
cur_content = "<br><br><a href=#top>Go to the top</a><hr>\n"
f_sp.write("<br>"+cur_let+".&nbsp;\n")
for (speckey, specval) in list_specials:
  #curlink = "<br /><a href=\"#"+abivarname+"\">"+abivarname+"</a> : "+defi+"\n"
  while not speckey.lower().startswith(cur_let.lower()):
    cur_let = chr(ord(cur_let)+1)
    f_sp.write("<br>"+cur_let+".&nbsp;\n")
  curlink = "<a href=\"#"+speckey+"\">"+speckey+"</a>&nbsp;&nbsp;\n"
  f_sp.write(curlink)

  cur_content += "<br><font id=\"title\"><a name=\""+speckey+"\">"+speckey+"</a></font>\n"
  cur_content += "<br><font id=\"text\">\n"
  cur_content += "<p>\n"+doku2html(make_links(specval,speckey,list_all_vars,list_chars,cur_specials))+"\n"
  cur_content += "</font>"
  cur_content += "<br><br><a href=#top>Go to the top</a><hr>\n"

f_sp.write(cur_content)

if make_topics_visible==0:
  cur_content += "<script type=\"text/javascript\" src=\""+js_path+"list_internal_links.js\"> </script>\n\n"
else:
  cur_content += "<script type=\"text/javascript\" src=\""+js_path+"list_internal_links_incl_topics.js\"> </script>\n\n"
cur_content += "</body>\n"
cur_content += "</html>"

f_sp.write(cur_content)
f_sp.close()
print("File specials.html has been written ...")

# Allvariables file : complete the content, then write the file and close it
if make_topics_visible==0:
  output += "<script type=\"text/javascript\" src=\""+js_path+"list_internal_links.js\"> </script>\n\n"
else:
  output += "<script type=\"text/javascript\" src=\""+js_path+"list_internal_links_incl_topics.js\"> </script>\n\n"
output += "</body>\n"
output += "</html>"

file_html = 'html_automatically_generated/allvariables.html'
f_html = open(file_html,'w')
f_html.write(output)
f_html.close()
print("File allvariables.html has been written ...")

################################################################################
################################################################################

# Generate the different "topic" files

################################################################################
# Constitute the section 3 "Related input variables" for all topic files. 
# This sec3 is stored, for each topic_name, in topic_sec3[topic_name]

topic_sec3 = dict()
topic_class_sec3 = dict()
found = dict()
foundvar = dict()
abivarname_dic = dict()
abivarname2_dic = dict()
section_dic = dict()
mnemonics_dic = dict()

for i, topic in enumerate(topics):
  topic_name = topic.topic_name

  topic_class_sec3[topic_name] = ""
  topic_sec3[topic_name] = ""
  found[topic_name] = 0
  foundvar[topic_name] = 0
  abivarname_dic[topic_name] = ""
  abivarname2_dic[topic_name] = ""
  section_dic[topic_name] = ""
  mnemonics_dic[topic_name] = ""

for (tclasskey, tclassval) in list_topics_class:

  for topic_name, value in topic_class_sec3.items():
    topic_class_sec3[topic_name] = "<p>"+tclassval+"<p>"

  if debug == 1:
    print("\nWork on "+tclasskey+"\n")

  for i, var in enumerate(variables):
    foundvar[topic_name] = 0
    try:
      if var.topics is not None:
        #print("\nabivarname:",var.abivarname)  
        #print("topics:",var.topics)  
        vartopics=var.topics
        topics_name_class = vartopics.split(',')
        for i, topic_name_class in enumerate(topics_name_class):
          name_class = topic_name_class.split('_')
          #print("tclasskey,name_class[1].strip()",tclasskey,name_class[1].strip()) 
          if tclasskey==name_class[1].strip() :
            topic_name=name_class[0].strip()
            found[topic_name] = 1
            foundvar[topic_name] = 1
            abivarname_dic[topic_name]=var.abivarname
            abivarname2_dic[topic_name]=var.abivarname
            if var.characteristics is not None and '[[INTERNAL_ONLY]]' in var.characteristics:
              abivarname2_dic[topic_name] = '%'+var.abivarname
            section_dic[topic_name]=var.section
            mnemonics_dic[topic_name]=var.mnemonics
    except:
      if debug==1 :
       print(" No topics for abivarname "+var.abivarname) 

    for i, topic in enumerate(topics):
      topic_name=topic.topic_name

      #print("topic_name: "+topic_name,"foundvar: ",foundvar[topic_name])

      if foundvar[topic_name] == 1:
        # Constitute the line of information related to one input variable
        foundvar[topic_name]=0
        abivarname=abivarname_dic[topic_name]
        abivarname2=abivarname2_dic[topic_name]
        section=section_dic[topic_name]
        mnemonics=mnemonics_dic[topic_name]
        topic_class_sec3[topic_name] += "... <a href=\""+section+".html#"+abivarname+"\">"+abivarname2+"</a>   "
        topic_class_sec3[topic_name] += "["+mnemonics+"]<br>\n"

        if debug==1 :
          print("topic_name:"+topic_name)
          print("topic_class_sec3[topic_name]:")
          print(topic_class_sec3[topic_name])

  for i, topic in enumerate(topics):
    topic_name=topic.topic_name
    if found[topic_name] == 1:
      found[topic_name]=0
      topic_sec3[topic_name] = topic_sec3[topic_name] + topic_class_sec3[topic_name]

################################################################################
# Constitute the section 4 "Selected input files" for all topic files.
# This sec4 is stored, for each topic_name, in topic_sec4[topic_name]

topic_sec4 = dict()
topic_class_sec4 = dict()
found = dict()
for i, var in enumerate(topics):
  topic_name = var.topic_name
  topic_class_sec4[topic_name] = ""
  topic_sec4[topic_name] = ""
  found[topic_name] = 0
 
topic_header_sec4="\n\n<p>&nbsp; \n<HR ALIGN=left> \n<p> <a name=\"4\">&nbsp;</a>\n<h3><b> 4. Selected input files.</b></h3>\n\n\n"
topic_header_sec4+="The user can find some related example input files in the ABINIT package in the directory /tests, or on the Web:\n"

# Create a dictionary to contain the list of tests for each topics
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
  found[topic_name] = 1
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
    topic_class_sec4[topic_name]+= line

for topic_name, value in found.items():
  if found[topic_name] == 1:
    topic_sec4[topic_name] = topic_header_sec4 + topic_class_sec4[topic_name]
    topic_sec4[topic_name] += "<br>\n"


################################################################################
# Assemble the "topic" files : secs 1, 2, 3 and then possibly 4 and 5

# For each "topic" file 
for topic_name, content in topic_sec3.items():
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

  f_topic.write(topic_header_varX)

# Write Sec. 3 
  f_topic.write(topic_sec3[topic_name])

# Write Sec. 4
  f_topic.write(topic_sec4[topic_name])

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

################################################################################
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
