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
from doc.pymods.lib_to_assemble_html import *

debug = 0

# Relative path from HTML files
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
################################################################################

# Parsing section, also preliminary treatment of list of topics in tests and bibliography 

################################################################################
 

bibyml=read_yaml_file(bib_ori+"/bibhtml.yml")
helps=read_yaml_file(help_ori+"/helps.yml")
lessons=read_yaml_file(tuto_ori+"/lessons.yml")
list_of_topics=read_yaml_file(topics_ori+"/list_of_topics.yml")
varfiles=read_yaml_file(invars_ori+"/varfiles.yml")
tests_dirs=read_yaml_file(topics_ori+"/tests_dirs.yml")
theorydocs=read_yaml_file(theory_ori+"/theorydocs.yml")
variables=read_yaml_file(invars_ori+"/abinit_vars.yml")

with open(bib_ori+'/abiref.bib')  as bibtex_file:
  bibtex_str = bibtex_file.read()

################################################################################
# Parse the ABINIT input files, in order to find the possible topics to which they are linked -> topics_in_tests

try :
  rm_cmd = "rm "+topics_gen+"/topics_in_tests.yml"
  retcode = os.system(rm_cmd)
except :
  if debug==1 :
    print("rm "+topics_gen+"/topics_in_tests.yml failed")
    print("the file was likely non existent")

for tests_dir in tests_dirs :
  grep_cmd = "grep topics tests/%s/Input/*.in > %s/topics_in_tests.txt"%(tests_dir,topics_gen)
  retcode = os.system(grep_cmd)
  if retcode == 0 :
    sed_cmd = "sed -e 's/^/- /' %s/topics_in_tests.txt >> %s/topics_in_tests.yml"%(topics_gen,topics_gen)
    retcode = os.system(sed_cmd)

topics_in_tests=read_yaml_file(topics_gen+"/topics_in_tests.yml")
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
# Collect the link seeds

allowed_link_seeds={}

# Groups of seeds
for var in variables:
  allowed_link_seeds[var.abivarname]="input_variable in "+var.varfile
for item in list_chars:
  allowed_link_seeds[item]="characteristic"
for (specialkey,specialval) in list_specials:
  allowed_link_seeds[specialkey]="special"
for i, varfile_info in enumerate(varfiles):
  varfile = varfile_info.name
  allowed_link_seeds[varfile]="varfile"
for i, lesson_info in enumerate(lessons):
  lesson = lesson_info.name
  allowed_link_seeds[lesson]="lesson"
for i, theory_info in enumerate(theorydocs):
  theorydoc = theory_info.name
  allowed_link_seeds[theorydoc]="theorydoc"
for i, help_info in enumerate(helps):
  helpfile = help_info.name
  allowed_link_seeds[helpfile]="helpfile"
for ref in bibtex_dics:
  ID=ref["ID"]
  allowed_link_seeds[ID]="bibID"

# Specific allowed seeds
allowed_link_seeds["allvariables"]="allvariables"

################################################################################
################################################################################

# The information needed to create automatically the links has been collected

################################################################################
# Initialization to constitute the body of the specials and var*.html files

all_vars = dict()
all_contents = dict()
for i, varfile_info in enumerate(varfiles):
  varfile = varfile_info.name
  all_vars[varfile] = []
  all_contents[varfile]= "<br><br><br><br><hr>\n"

################################################################################
# Constitute the body of information for the special parameters, stored for the appropriate varfile in all_contents[varfile]

cur_specials = []
for (specialkey,specialval) in list_specials:
  cur_specials.append(specialkey)

for (speckey, specval) in list_specials:
  backlink= ' &nbsp; <a href="../../%s/specials.html#%s">%s</a> &nbsp; ' %(invars_gen,speckey,speckey)
  cur_content = "<br><font id=\"title\"><a name=\""+speckey+"\">"+speckey+"</a></font>\n"
  cur_content += "<br><font id=\"text\">\n"
  cur_content += "<p>\n"+doku2html(make_links(specval,speckey,allowed_link_seeds,backlinks,backlink))+"\n"
  cur_content += "</font>"
  cur_content += "<br><br><a href=#top>Go to the top</a>\n"
  cur_content += "<B> | </B><a href=\"allvariables.html#top\">Complete list of input variables</a><hr>\n"
  #
  all_contents["specials"] = all_contents["specials"] + cur_content + "\n\n"

################################################################################
# Constitute the body of information for all variables, stored for the appropriate varfile in all_contents[varfile]

for i, var in enumerate(variables):
  # Constitute the body of information related to one input variable
  varfile = var.varfile
  all_vars[varfile].append([var.abivarname,var.mnemonics])
  cur_content = ""
  backlink=' &nbsp; <a href="../../%s/%s.html#%s">%s</a> &nbsp; ' %(invars_gen,varfile,var.abivarname,var.abivarname)

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
      cur_content += "<br><font id=\"characteristic\">Characteristic: "+make_links(chars,var.abivarname,allowed_link_seeds,backlinks,backlink)+"</font>\n"
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
          cur_content += '<a href="../../'+topics_gen+'/topic_'+name_class[0].strip()+'.html">'+name_class[0].strip()+'</a> '
        cur_content += "</font>\n"
    except:
      if debug==1 :
        print(" No topic_class for abivarname "+var.abivarname)
    # Variable type, including dimensions
    cur_content += "<br><font id=\"vartype\">Variable type: "+var.vartype
    if var.dimensions is not None:
      cur_content += make_links(format_dimensions(var.dimensions),var.abivarname,allowed_link_seeds,backlinks,backlink)
    if var.commentdims is not None and var.commentdims != "":
      cur_content += " (Comment: "+make_links(var.commentdims,var.abivarname,allowed_link_seeds,backlinks,backlink)+")"
    cur_content += "</font>\n" 
    # Default
    cur_content += "<br><font id=\"default\">"+make_links(format_default(var.defaultval),var.abivarname,allowed_link_seeds,backlinks,backlink)
    if var.commentdefault is not None and var.commentdefault != "":
      cur_content += " (Comment: "+make_links(var.commentdefault,var.abivarname,allowed_link_seeds,backlinks,backlink)+")"
    cur_content += "</font>\n" 
    # Requires
    if var.requires is not None and var.requires != "":
      cur_content += "<br><br><font id=\"requires\">\nOnly relevant if "+doku2html(make_links(var.requires,var.abivarname,allowed_link_seeds,backlinks,backlink))+"\n</font>\n"
    # Excludes
    if var.excludes is not None and var.excludes != "":
      cur_content += "<br><br><font id=\"excludes\">\nThe use of this variable forbids the use of "+doku2html(make_links(var.excludes,var.abivarname,allowed_link_seeds,backlinks,backlink))+"\n</font>\n"
    # Text
    cur_content += "<br><font id=\"text\">\n"
    cur_content += "<p>\n"+doku2html(make_links(var.text,var.abivarname,allowed_link_seeds,backlinks,backlink))+"\n"
    # End the section for one variable
    cur_content += "</font>\n\n"
    cur_content += "<br><br><a href=#top>Go to the top</a>\n"
    cur_content += "<B> | </B><a href=\"allvariables.html#top\">Complete list of input variables</a><hr>\n"
    #
    all_contents[varfile] = all_contents[varfile] + cur_content + "\n\n"
  except AttributeError as e:
    print(e)
    print('For variable : ',abivarname)

################################################################################
# Generate the files that document all the variables (all such files : var* as well as all and special).

suppl_components={}

# Store the default informations
for i, varfile_info in enumerate(varfiles):
  if varfile_info.name.strip()=="default":
    varfile_info_default=varfile_info

# Generate each "normal" varfile file : build the missing information (table of content), assemble the content, apply global transformations, then write.
for i, varfile_info in enumerate(varfiles):
  varfile = varfile_info.name
  if varfile=="default":
    continue

  #Generate the body of the table of content 
  cur_let = 'A'
  toc_body = " <br>"+cur_let+".\n"
  if varfile=="allvariables":
    for i, var in enumerate(variables):
      while not var.abivarname.startswith(cur_let.lower()):
        cur_let = chr(ord(cur_let)+1)
        toc_body += " <p>"+cur_let+".\n"
      abivarname=var.abivarname
      if var.characteristics is not None and '[[INTERNAL_ONLY]]' in var.characteristics:
        abivarname = '%'+abivarname
      curlink = " <a href=\""+var.varfile+".html#"+var.abivarname+"\">"+abivarname+"</a>&nbsp;&nbsp;\n"
      toc_body += curlink
  elif varfile == "specials":
    for (speckey, specval) in list_specials:
      while not speckey.lower().startswith(cur_let.lower()):
        cur_let = chr(ord(cur_let)+1)
        toc_body += " <br>"+cur_let+".\n"
      curlink = " <a href=\"#"+speckey+"\">"+speckey+"</a>&nbsp;&nbsp;\n"
      toc_body += curlink
  else:
    for abivarname,defi in all_vars[varfile]:
      while not abivarname.startswith(cur_let.lower()):
        cur_let = chr(ord(cur_let)+1)
        toc_body += " <br>"+cur_let+".\n"
      curlink = " <a href=\"#"+abivarname+"\">"+abivarname+"</a>&nbsp;&nbsp;\n"
      toc_body += curlink
  toc_body += "\n"

  suppl={"toc":toc_body , "content":all_contents[varfile]}
  suppl_components[varfile]=suppl

rc=assemble_html(varfiles,suppl_components,"input_variables","",allowed_link_seeds,backlinks)

################################################################################
################################################################################

# Generate the different "topic" files

################################################################################
# Constitute the component "Related input variables" for all topic files. 
# This component in input variables is stored, for each topic_name, in topic_invars[topic_name]

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
            topic_invars[topic_name] += '... <a href="../../'+invars_gen+'/'+var.varfile+'.html#'+var.abivarname+'">'+abivarname+'</a>   '
            topic_invars[topic_name] += "["+var.mnemonics+"]<br>\n"
    except:
      if debug==1 :
       print(" No topics for abivarname "+var.abivarname) 

################################################################################
# Constitute the component "Selected input files" for all topic files.
# This component on input files in topics  is stored, for each topic_name, in topic_infiles[topic_name]

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
      dirname=file_split[1] 
      testname=file_split[3] 
      if dirname not in dir.keys():
        dir[dirname] = []
      dir[dirname].append(testname)
    for dirname, testnames in dir.items():
      line="<p> tests/"+dirname+"/Input: "
      for testname in testnames:
        line+="<a href=\"../../tests/"+dirname+"/Input/"+testname+"\">"+testname+"</a> \n"
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
default_topic_yml=read_yaml_file(topics_ori+"/default_topic.yml")
default_topic=default_topic_yml[0]

# For each "topic" file
for topic_name in list_of_topics:
  path_yml=topics_ori+"/topic_"+topic_name+".yml"
  print("Will use file "+path_yml+" to initiate the topic "+topic_name+" ... ",end="")
  topic_yml=read_yaml_file(path_yml)
  topic=topic_yml[0]

  #Mention it in the table of content of the file all_topics.html
  while not (topic_name.startswith(cur_let_all.lower()) or topic_name.startswith(cur_let_all.upper())):
    cur_let_all = chr(ord(cur_let_all)+1)
    toc_all = toc_all + "<p>"+cur_let_all+".\n"
  toc_all = toc_all + "<br><a href=\"topic_"+ topic_name + ".html\">" + topic_name + "</a> [How to "+topic.howto+"] &nbsp;&nbsp;\n"

  #Find the bibliographical references
  reflist=[]
  for j in ["introduction","examples"] :
    try :
      extract_j=getattr(topic,j).strip()
    except :
      extract_j=""
    linklist=re.findall("\\[\\[([a-zA-Z0-9_ */<>]*)\\]\\]",extract_j,flags=0)
    for ref in linklist:
      m=re.search("\d{4}",ref,flags=0)
      if m!=None:
        reflist.append(ref)

  reflist=list(set(reflist))
  reflist.sort()

  topic_refs=""
  for (i,ID) in enumerate(reflist):
    topic_refs+="<br> [["+ID+"]] "+reference_dic[ID]+"<br> \n"
  topic_refs+="<p>"
  topic_refs=bibtex2html(topic_refs)

  #Generate the table of content
  item_toc=0
  item_list=[]
  title={ "introduction":"Introduction." , "examples":"Example(s)", "tutorials":"Related tutorials." , 
          "input_variables":"Related input variables." , "input_files":"Selected input files." , "references":"References."}
  sec_number={}
  toc=" <h3><b>Table of content: </b></h3> \n <ul> "
  for j in ["introduction","examples","tutorials","input_variables","input_files","references"] :
    sec_number[j]="0"
    try :
      extract_j=getattr(topic,j).strip()
    except :
      extract_j=""
    if (extract_j != "" and extract_j!= "default") or (j=="input_variables" and topic_invars[topic_name]!="") or (j=="input_files" and topic_infiles[topic_name]!="") or (j=="references" and reflist!=[]):
      item_toc += 1
      item_num="%d" % item_toc
      sec_number[j]=item_num
      toc += '<li><a href="topic_'+topic_name+'.html#'+item_num+'">'+item_num+'</a>. '+title[j]
      
  toc+= "</ul>"

  #Generate a first version of the html file, in the order "header" ... up to the "end"
  #Take the info from the component "default" if there is no information on the specific component provided in the yml file.
  topic_html=""
  for j in ["header","title","subtitle","copyright","links","toc","introduction","examples","tutorials","input_variables","input_files","references","links","end"]:
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
    elif j == "references":
      if sec_number[j]!="0" :
        topic_html+= '\n&nbsp; \n<HR ALIGN=left> \n<a name=\"'+sec_number[j]+'\">&nbsp;</a>\n<h3><b>'+sec_number[j]+'. '+title[j]+'</b></h3>\n\n\n'
        topic_html+= topic_refs
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

  dir_root="topics"
  name_root="topic"
  rc=finalize_html(topic_html,topic,dir_root,name_root,allowed_link_seeds,backlinks)

################################################################################
# Generate the file with the list of names of different "topic" files

#Generate a first version of the html file, in the order "header" ... up to the "end"
#Take the info from the component "default" if there is no information on the specific component provided in the yml file.
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

dir_root="topics"
name_root=""
rc=finalize_html(all_topics_html,default_topic,dir_root,name_root,allowed_link_seeds,backlinks)

################################################################################
################################################################################

# Automatic translation
# of the *_help.html files to help_*.yml files

################################################################################

activate_translation=1
if activate_translation==1:

  docs=read_yaml_file(help_ori+"/helps.yml")

  for i, doc_info in enumerate(docs):

    # Skip the default component
    name = doc_info.name
    if name=="default":
      break
    path_doc_html="users/"+name+"_help.html"
    path_doc_yml=help_ori+"/help_"+name+".yml"
    print("Will use file "+path_doc_html+" to build the file "+path_doc_yml+" ... ",end="")

    f_doc_html=open(path_doc_html,"r")
    doc_html=f_doc_html.readlines()

    doc_yml=""
    doc_yml+="# This YAML file contains the introduction as well as the body (including the table of content) of the html help document.\n"
    doc_yml+="# In order to modify the other parts, modify the file helps.yml .\n"
    doc_yml+="# This is the introduction ...\n"
    doc_yml+="intro : |\n"

    body_header=""
    body_header+="# This is the body, including the table of content ...\n"
    body_header+="body : |\n"

    intro=0
    body=0
    for line in doc_html:
      if "<!--" in line and "-->" in line:
        if "begin" in line :
          if intro==1 or body==1:
            raise ValueError("(intro,body)=(%s,%s)"%(intro,body))
          if "intro" in line :
            intro=1
            continue
          if "body" in line:
            body=1
            doc_yml+=body_header
            continue
        if "end" in line and "intro" in line:
          if intro==0 or body==1:
            raise ValueError("(intro,body)=(%s,%s)"%(intro,body))
          intro=0
        if "end" in line and "body" in line:
          if intro==1 or body==0:
            raise ValueError("(intro,body)=(%s,%s)"%(intro,body))
          body=0

      if intro+body==1 :
        #The line must be written, but must possibly perform changes:
        if "<a href=" in line:
          # Stabilize the own reference
          string_old='href="'+name+'_help.html'
          string_new='href="'
          line=line.replace(string_old,string_new)
          string_old='href="./'+name+'_help.html'
          string_new='href="'
          line=line.replace(string_old,string_new)
          # Correct the references to the other files in the tutorial directory (transient measure in case of the "lesson_" files)
          string_old='href="lesson_'
          string_new='href="../../tutorial/lesson_'
          line=line.replace(string_old,string_new)
          string_old='href="./lesson_'
          line=line.replace(string_old,string_new)
          #string_old='href="theory_'
          #string_new='href="../theory_'
          #line=line.replace(string_old,string_new)
          #string_old='href="./theory_'
          #line=line.replace(string_old,string_new)
          #string_old='href="welcome'
          #string_new='href="../welcome'
          #line=line.replace(string_old,string_new)
          #string_old='href="./welcome'
          #line=line.replace(string_old,string_new)
          # Create automatically the new links for the input variables
          if "html_automatically_generated" in line:
            # See whether one variable appear
            for i, var in enumerate(variables):
              if var.abivarname in line:
                varname = var.abivarname
                varfile = var.varfile
                string_old='<a href="../input_variables/html_automatically_generated/%s.html#%s" target="kwimg">%s</a>'%(varfile,varname,varname)
                string_new="[["+varname+"]]"
                line=line.replace('"'+string_old+'"',string_new)
                line=line.replace(string_old,string_new)
                # Slight variation
                string_old='<a href="../input_variables/html_automatically_generated/%s.html#%s" target="kwimg">%s</a>'%(varfile,varname,varname)
                line=line.replace('"'+string_old+'"',string_new)
                line=line.replace(string_old,string_new)
            # Otherwise, correct the path
            string_old='href="../input_variables/html_automatically_generated'
            string_new='href="../../input_variables/generated_files'
            line=line.replace(string_old,string_new)
          if "users" in line:
            string_old='href="../users/'
            string_new='href="../../users/'
            line=line.replace(string_old,string_new)

        #string_old='src="theory'
        #string_new='src="../documents/theory'
        #line=line.replace(string_old,string_new)
        #string_old='src="./theory'
        #line=line.replace(string_old,string_new)
        #string_old='src=./theory'
        #string_new='src=../documents/theory'
        #line=line.replace(string_old,string_new)

        doc_yml+="  "+line

    #print("")
    #print(" doc_yml :",path_doc_yml)
    #index=36700
    #print(" doc_yml[%s,%s] :"%(index,index+100))
    #print(doc_yml[index:index+100])
    #print("")

    # Write the finalized html file
    f_doc_yml=open(path_doc_yml,"w")
    f_doc_yml.write(doc_yml)
    f_doc_yml.close()
    print("File %s written ..."%path_doc_yml )

################################################################################
################################################################################

# Assemble the html files to be generated from the yml information.
# In order : tutorial, files lessons_*
#            theory, files theory_*
#            users,  files *_help
 
################################################################################

suppl_components={}
returncode=assemble_html(lessons,suppl_components,"tutorial","lesson",allowed_link_seeds,backlinks)
returncode=assemble_html(theorydocs,suppl_components,"theory","theory",allowed_link_seeds,backlinks)
returncode=assemble_html(helps,suppl_components,"users","help",allowed_link_seeds,backlinks)

################################################################################
################################################################################

# Come back to the bibliography

################################################################################
# Treat the links within the "introduction" of the acknowledgment component first.

backlink= ' &nbsp; <a href="../../%s/acknowledgments.html">acknowledgments.html</a> &nbsp; ' %(bib_gen)
for i, bibyml_info in enumerate(bibyml):
  if bibyml_info.name.strip()=="acknowledgments":
    bibyml_intro=bibyml_info.introduction
    bibyml_ack_intro = doku2html(make_links(bibyml_intro,None,allowed_link_seeds,backlinks,backlink))

################################################################################
# Write an ordered bib file, that allows to update the original one.
# Constitute also the content of the ABINIT bibtex and bibliography files.

bibtex_content=""
bibliography_content=""
lines_txt=""
cur_let = 'A'
alphalinks="\n \n <hr> Go to "
for i in string.ascii_uppercase:
  alphalinks+=('<a href=#%s>%s</a> ')%(i,i)
alphalinks+="\n \n"
bibliography_content+=('<a id="%s"></a>')%(cur_let)+alphalinks+('<hr><hr><h2>%s</h2> \n \n')%(cur_let)
for ref in bibtex_dics:
  entrytype=ref["ENTRYTYPE"]
  ID=ref["ID"]
  list_backlinksID=backlinks[ID].split(";;")
  for (i,link) in enumerate(list_backlinksID):
     list_backlinksID[i]=link.strip()
  backlinksID=set(list_backlinksID)
  line=("@%s{%s,%s") %(entrytype,ID,ref['body'])
  lines_txt+= line
  bibtex_content+= ('<hr><a id="%s">%s</a> \n <pre>' ) %(ID,ID)
  bibtex_content+= line+'</pre> \n'
  if ID[0]>cur_let:
    cur_let=ID[0]
    bibliography_content+=('<a id="%s"></a>')%(cur_let)+alphalinks+('<hr><hr><h2>%s</h2> \n \n')%(cur_let)
  bibliography_content+= ('<hr><a id="%s">[%s]</a> (<a href="./bibtex.html#%s">bibtex</a>)\n <br> %s \n') %(ID,ID,ID,reference_dic[ID])
  nlink=0
  for link in backlinksID:
    if len(link)!=0:
      if nlink==0: 
        bibliography_content+= "<br> Referred to in " 
        nlink=1
      bibliography_content+= link

# Open, write and close the txt file
file_txt = bib_gen+'/ordered_abiref.bib'
f_txt = open(file_txt,'w')
f_txt.write(lines_txt)
f_txt.close()
print("File %s has been written ..." %file_txt)

################################################################################
#Global operation on the bibliography html file : conversion from bibtex notation to html notation.
#This should cover most of the cases. 
 
bibliography_content=bibtex2html(bibliography_content)

################################################################################
# Assemble the html files in the bibliography directory

suppl={"introduction":bibyml_ack_intro}
suppl_components={"acknowledgments":suppl}
suppl={"content":bibtex_content}
suppl_components['bibtex']=suppl
suppl={"content":bibliography_content}
suppl_components['bibliography']=suppl

rc=assemble_html(bibyml,suppl_components,"bibliography","",allowed_link_seeds,backlinks)

################################################################################

print("Work done !")

################################################################################
