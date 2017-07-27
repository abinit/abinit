#! /usr/bin/env python 
"""
   This script generates most of the ABINIT documentation for the ABINIT Web site, in HTML.

   Usage : python generate_doc.py [-h|--help]
   (No options are allowed at present, except the help ones)

   This script :
   (1) reads information from (i) a set of YAML files, (ii) a bibtex file, 
        (iii) input files contained in tests/*/Input 
       (files (i) and (ii) are contained in subdirectories */origin_files);
   (2) performs some checks; 
   (3) establishes intermediate dictionaries and databases;
   (4) expands special strings, of the form "[[tag]]", to create HTML links;
   (5) creates the needed HMTL files (all contained in subdirectories */generated_files).

   The result of the expansion of "[[tag]]" will depend on the tag. If valid, the URL will of course
   point to the correct location. The string that will be seen on screen will usually
   be the string stripped from the pairs of square brackets, but there is an exception, see later.
   Also, only selected classes of tags can be expanded. Other strings will be declared "FAKE_LINK" on the screen.
   
   A. If tag is an ABINIT input variable, mentioned in the database ~abinit/doc/input_variables/origin_files/abinit_vars.yml,
      or equivalently in the list http://www.abinit.org/inpvars,
      "[[tag]]" will appear on screen as "tag" without the two pairs of squarebrackets.
      Example: "[[acell]]" will apear on screen as "acell".

   B. If tag is a bibliographical reference, mentioned in the database ~abinit/doc/bibliography/origin_files/abiref.bib,
      "[[tag]]" will appear on screen as "[tag]" with onoly one of the two pairs of squarebrackets.
      Note that such tag must have the form of the name of the first authors, with the first letter uppercase and the other lower case,
      followed by the year of publication (four digits) and possibly a letter if more than one article would have the same tag.
      Example: "[[Kohn1965]]" will appear on screen as "[Kohn1965]".
               "[[Amadon2008b]]" will appear on screen as "[Amadon2008b]".

   C. The tags that starts with "lesson_", "topic_", "theorydoc_", "var", "allvar", "help_" and corresponds to one of the existing
      lessons of the tutorial, or one of the topics, theorydocs, input variables files, or help files, are allowed tags. 
      Most of them will appear on screen as "tag" without the two pairs of squarebrackets, EXCEPT the
      "help_XYZ" ones, that will appear as "XYZ help file".
      Examples: "[[lesson_base1]]" will appear on screen as "lesson_base1"
                "[[help_new_user]]" will appear on screen as "new_user help file"
"""

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

################################################################################
################################################################################

# Generic information


################################################################################

cmdline_params = sys.argv[1:]
if cmdline_params != [] :
  if cmdline_params[0] == "-h" or cmdline_params[0] == "--help" :
    print(__doc__)
  else :
    print (" Unrecognized option. Abort.")
  sys.exit()

################################################################################
################################################################################

# Parsing section, also preliminary treatment of list of topics in tests and bibliography 

################################################################################

list_infos_dir=[]
list_infos_dir.append({"dir_name":"bibliography","root_filname":"",
                                                    "yml_files":["bibfiles"]})
list_infos_dir.append({"dir_name":"input_variables","root_filname":"",
                                                    "yml_files":["abinit_vars","characteristics","list_specials","varfiles"]})
list_infos_dir.append({"dir_name":"theory","root_filname":"theorydoc",
                                                    "yml_files":["theorydocs"]})
list_infos_dir.append({"dir_name":"topics","root_filname":"topic",
                                                    "yml_files":["default_topic","list_of_topics","list_tribes","tests_dirs"]})
list_infos_dir.append({"dir_name":"tutorial","root_filname":"lesson",
                                                    "yml_files":["lessons"]})
list_infos_dir.append({"dir_name":"users","root_filname":"help",
                                                    "yml_files":["helps"]})
msgs={"bibfiles"       :"as database input file for the list of generated files in the bibliography directory ...",
      "abinit_vars"    :"as database input file for the input variables and their characteristics ...",
      "characteristics":"as database input file for the list of allowed characteristics ...",
      "list_specials"  :"as database input file for the list of allowed special keywords ...",
      "varfiles"       :"as database input file for the list of varfiles ...",
      "theorydocs"     :"as database input file for the list of theory documents ...",
      "default_topic"  :"to initialize the topic html files with default values ...",
      "list_of_topics" :"as database input file for the list of topics ...",
      "list_tribes"    :"as database input file for the list of tribes ...",
      "tests_dirs"     :"as database input file for the list of directories in which automatic test input files are present ...",
      "lessons"        :"as database input file for the list of lessons ...",
      "helps"          :"as database input file for the list of help files in the users directory ..."}
      
path_file='bibliography/origin_files/abiref.bib'
with open(path_file) as f:
  print("Read "+path_file+" as database input file for the bibliography references ...")
  bibtex_str = f.read()

yml_in={}
for infos_dir in list_infos_dir:
  yml_files=infos_dir["yml_files"]
  for yml_file in yml_files:
    path_ymlfile="%s/origin_files/%s.yml"%(infos_dir["dir_name"],yml_file)
    print("Read "+path_ymlfile+" "+msgs[yml_file])
    yml_in[yml_file] =read_yaml(path_ymlfile)

# These ones are quite often used, so copy them ...
abinit_vars=yml_in["abinit_vars"]
list_specials=yml_in["list_specials"]
varfiles=yml_in["varfiles"]
list_of_topics=yml_in["list_of_topics"]
  
################################################################################
# Parse the ABINIT input files, in order to find the possible topics to which they are linked -> topics_in_tests
# Also constitute the list of allowed links to tests files.

try :
  rm_cmd = "rm topics/generated_files/topics_in_tests.txt"
  retcode = os.system(rm_cmd)
except :
  if debug==1 :
    print(rm_cmd+"failed")
    print("the file was likely non existent")
try :
  rm_cmd = "rm topics/generated_files/topics_in_tests.yml"
  retcode = os.system(rm_cmd)
except :
  if debug==1 :
    print(rm_cmd+"failed")
    print("the file was likely non existent")

allowed_links_in_tests=[]
for tests_dir in yml_in["tests_dirs"] :

  grep_cmd = "grep topics tests/%s/Input/*.in > topics/generated_files/topics_in_tests.txt"%(tests_dir)
  retcode = os.system(grep_cmd)
  if retcode == 0 :
    sed_cmd = "sed -e 's/^/- /' topics/generated_files/topics_in_tests.txt >> topics/generated_files/topics_in_tests.yml"
    retcode = os.system(sed_cmd)

  # Allowed links
  path_dir_input="tests/%s/Input"%(tests_dir)
  list_files=os.listdir(path_dir_input)
  for file in list_files:
    allowed_links_in_tests.append(path_dir_input+'/'+file)
  if tests_dir[:4] == 'tuto':
    path_dir_refs="tests/%s/Refs"%(tests_dir)
    list_files=os.listdir(path_dir_refs)
    for file in list_files:
      allowed_links_in_tests.append(path_dir_refs+'/'+file)

path_ymlfile="topics/generated_files/topics_in_tests.yml"
print("Generated "+path_ymlfile+", to contain the list of automatic test input files relevant for each topic ...")
topics_in_tests=read_yaml(path_ymlfile)
try :
  rm_cmd = "rm topics/generated_files/topics_in_tests.txt"
  retcode = os.system(rm_cmd)
except :
  if debug==1 :
    print(rm_cmd+"failed")
    print("the file was likely non existent")

if debug==1 :
  print(" topics_in_tests :")
  print(topics_in_tests)

################################################################################
# Constitutes a list of bib items, each being a dictionary

bibtex_dics=[]
bibtex_items=bibtex_str.split('@')
empty=bibtex_items.pop(0)
listIDs=[]

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
  if not entrytype in ["article","book","incollection","phdthesis","mastersthesis","misc"]:
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
  listIDs.append(id_upper_lower)

  #Store the remaining, for later possible reordering, without any of the later treatments.
  item_dic['body']=item[1]

  item[1]=item[1].replace('optdoi','doi')
  item[1]=item[1].replace('opturl','url')
  item[1]=item[1].replace('optURI','url')
  item[1]=item[1].replace('adsurl','url')

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
  note=""
  howpublished=""
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
    elif key=='note':
      note=ref.values()[position].strip()
    elif key=='howpublished':
      howpublished=ref.values()[position].strip()
    position+=1

  # Reformat the list of authors, starting with the initials.
  author=reformat_namelist(author)
  editor=reformat_namelist(editor)

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
  elif ENTRYTYPE=="mastersthesis":
    formatted=' %s "%s", Masters thesis (%s, %s, %s).' % (author,title,school,address,year)
  elif ENTRYTYPE=="book":
    if volume!="":
      formatted=' %s "%s", vol.%s (%s, %s, %s).' % (author,title,volume,publisher,address,year)
    else:
      formatted=' %s "%s" (%s, %s, %s).' % (author,title,publisher,address,year)
  elif ENTRYTYPE=="misc":
    formatted=' %s "%s", %s, %s (%s).' % (author,title,note,howpublished,year)

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
# Check that there is no conflict of references, then order the references
listIDs_sorted=sorted(listIDs)
for (i,item) in enumerate(listIDs):
  if i==0:
    continue
  if listIDs[i]==listIDs[i-1]:
    print(" Detected two bibtex items with the same ID:")
    print(" Item 1 has ID:",listIDs[i-1])
    print(" Item 2 has ID:",listIDs[i])
    print(" Please fix this issue ...")
    sys.exit()

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
file_txt = 'bibliography/generated_files/abiref.txt'
f_txt = open(file_txt,'w')
f_txt.write(lines_txt)
f_txt.close()
print("File %s written ..." %file_txt)

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
file_yml = 'bibliography/generated_files/abiref.yml'
f_yml = open(file_yml,'w')
f_yml.write(lines_yml)
f_yml.close()
print("File %s written ..." %file_yml)

################################################################################
# Collect the link seeds

allowed_link_seeds={}

# Groups of seeds
for var in abinit_vars:
  allowed_link_seeds[var.abivarname]="input_variable in "+var.varfile
for item in yml_in["characteristics"]:
  allowed_link_seeds[item]="characteristic"
for (specialkey,specialval) in list_specials:
  allowed_link_seeds[specialkey]="special"
for i, varfile_info in enumerate(varfiles):
  varfile = varfile_info.name
  allowed_link_seeds[varfile]="varfile"
for i, lesson_info in enumerate(yml_in["lessons"]):
  lesson = lesson_info.name
  allowed_link_seeds["lesson_"+lesson]="lesson"
for i, theory_info in enumerate(yml_in["theorydocs"]):
  theorydoc = theory_info.name
  allowed_link_seeds["theorydoc_"+theorydoc]="theorydoc"
for i, help_info in enumerate(yml_in["helps"]):
  help = help_info.name
  allowed_link_seeds["help_"+help]="help"
for ref in bibtex_dics:
  ID=ref["ID"]
  allowed_link_seeds[ID]="bibID"
for file in allowed_links_in_tests:
  allowed_link_seeds[file]="in_tests"

# Specific allowed seeds
allowed_link_seeds["varset_allvars"]="varset_allvars"

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
  backlink= ' &nbsp; <a href="../../input_variables/generated_files/specials.html#%s">%s</a> &nbsp; ' %(speckey,speckey)
  cur_content = "<br><font id=\"title\"><a name=\""+speckey+"\">"+speckey+"</a></font>\n"
  cur_content += "<br><font id=\"text\">\n"
  cur_content += "<p>\n"+make_links(specval,speckey,allowed_link_seeds,backlinks,backlink)+"\n"
  cur_content += "</font>"
  cur_content += "<br><br><a href=#top>Go to the top</a>\n"
  cur_content += "<B> | </B><a href=\"varset_allvars.html#top\">Complete list of input variables</a><hr>\n"
  #
  all_contents["specials"] = all_contents["specials"] + cur_content + "\n\n"

################################################################################
# Constitute the body of information for all variables, stored for the appropriate varfile in all_contents[varfile]

for i, var in enumerate(abinit_vars):
  # Constitute the body of information related to one input variable
  varfile = var.varfile
  all_vars[varfile].append([var.abivarname,var.mnemonics])
  cur_content = ""
  backlink=' &nbsp; <a href="../../input_variables/generated_files/%s.html#%s">%s</a> &nbsp; ' %(varfile,var.abivarname,var.abivarname)

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
        topics_name_tribe = vartopics.split(',')
        for i, topic_name_tribe in enumerate(topics_name_tribe):
          name_tribe = topic_name_tribe.split('_')
          cur_content += '<a href="../../topics/generated_files/topic_'+name_tribe[0].strip()+'.html">'+name_tribe[0].strip()+'</a> '
        cur_content += "</font>\n"
    except:
      if debug==1 :
        print(" No topic_tribe for abivarname "+var.abivarname)
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
      cur_content += "<br><br><font id=\"requires\">\nOnly relevant if "+make_links(var.requires,var.abivarname,allowed_link_seeds,backlinks,backlink)+"\n</font>\n"
    # Excludes
    if var.excludes is not None and var.excludes != "":
      cur_content += "<br><br><font id=\"excludes\">\nThe use of this variable forbids the use of "+make_links(var.excludes,var.abivarname,allowed_link_seeds,backlinks,backlink)+"\n</font>\n"
    # Text
    cur_content += "<br><font id=\"text\">\n"
    cur_content += "<p>\n"+make_links(var.text,var.abivarname,allowed_link_seeds,backlinks,backlink)+"\n"
    # End the section for one variable
    cur_content += "</font>\n\n"
    cur_content += "<br><br><a href=#top>Go to the top</a>\n"
    cur_content += "<B> | </B><a href=\"varset_allvars.html#top\">Complete list of input variables</a><hr>\n"
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

  scriptTab = "\n\
<input type=\"text\" id=\"InputSearch\" onkeyup=\"searchInput()\" onClick=\"searchInput()\" onblur=\"defaultClick()\" placeholder=\"Search\">\n\
"
  alphalinks="\n \n <div class=\"TabsLetter\">"
  for i in string.ascii_uppercase:
    alphalinks+=('<a class=\"TabLetterLink" href=\"#%s\" onClick=\"openLetter(event,\'%s\')\" id="click%s">%s</a> ')%(i,i,i,i)
  alphalinks+="</div>\n \n"

  #Generate the body of the table of content 
  cur_let = 'A'
  toc_body=""
  if varfile=="varset_allvars":
    toc_body+=scriptTab+alphalinks
  else :
    toc_body += " <br><a id='%s'></a>"%(cur_let)+cur_let+".\n"

  if varfile=="varset_allvars":
    toc_body += " <ul id=\"Letters\">\n"
    toc_body += " <li>\n<ul id=\"%s\" class=\"TabContentLetter\">\n"%(cur_let)
    toc_body += " <li class=\"HeaderLetter\">%s</li>\n"%(cur_let)
    for i, var in enumerate(abinit_vars):
      while not var.abivarname.startswith(cur_let.lower()):
        cur_let = chr(ord(cur_let)+1)
        toc_body += " </ul></li>\n<li>\n<ul id=\"%s\" class=\"TabContentLetter\">\n"%(cur_let)
        toc_body += " <li class=\"HeaderLetter col-s-6 col-m-3 col-l-2 col-xl-2 col-xxl-1\">%s</li>\n"%(cur_let)
      abivarname=var.abivarname
      if var.characteristics is not None and '[[INTERNAL_ONLY]]' in var.characteristics:
        abivarname = '%'+abivarname
      curlink = " <li class=\"col-s-6 col-m-3 col-l-2 col-xl-2 col-xxl-1\"><a href=\""+var.varfile+".html#"+var.abivarname+"\">"+abivarname+"</a></li>\n"
      toc_body += curlink
    toc_body += "</ul></li></ul>\n\
<script>\n\
defaultClick(true);\n\
</script>\n\
"
  elif varfile == "specials":
    for (speckey, specval) in list_specials:
      while not speckey.lower().startswith(cur_let.lower()):
        cur_let = chr(ord(cur_let)+1)
        toc_body += " <br><a id='%s'></a>"%(cur_let)+cur_let+".\n"
      curlink = " <a href=\"#"+speckey+"\">"+speckey+"</a>&nbsp;&nbsp;\n"
      toc_body += curlink
  else:
    for abivarname,defi in all_vars[varfile]:
      while not abivarname.startswith(cur_let.lower()):
        cur_let = chr(ord(cur_let)+1)
        toc_body += " <br><a id='%s'></a>"%(cur_let)+cur_let+".\n"
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

for (tribekey, tribeval) in yml_in["list_tribes"]:

  if debug == 1:
    print("\nWork on "+tribekey+"\n")

  for topic_name in list_of_topics:
    found[topic_name] = 0

  for i, var in enumerate(abinit_vars):
    try:
      if var.topics is not None:
        topics_name_tribe = var.topics.split(',')
        for i, topic_name_tribe in enumerate(topics_name_tribe):
          name_tribe = topic_name_tribe.split('_')
          if tribekey==name_tribe[1].strip() :
            topic_name=name_tribe[0].strip()
            if found[topic_name]==0 :
              topic_invars[topic_name] += "<p>"+tribeval+":<p>"
              found[topic_name] = 1
            abivarname=var.abivarname
            if var.characteristics is not None and '[[INTERNAL_ONLY]]' in var.characteristics:
              abivarname = '%'+abivarname
            topic_invars[topic_name] += '... <a href="../../input_variables/generated_files/'+var.varfile+'.html#'+var.abivarname+'">'+abivarname+'</a>   '
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
# Assemble the "topic" files 
# Also collect the keywords and howto 

default_topic=yml_in["default_topic"][0]
dic_keyword_name={}
dic_keyword_howto={}

# For each "topic" file
for topic_name in list_of_topics:
  path_ymlfile="topics/origin_files/topic_"+topic_name+".yml"
  print("Read "+path_ymlfile+" to initiate the topic '"+topic_name+"' ... ",end="")
  topic_yml=read_yaml(path_ymlfile) 
  topic=topic_yml[0]
  dic_keyword_name[topic.keyword]=topic_name
  dic_keyword_howto[topic.keyword]=topic.howto

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
  title={ "introduction":"Introduction." , "examples":"Example(s)", "tutorials":"Related lesson(s) of the tutorial." , 
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
# Note : they are ordered according to the keyword list

list_keywords=sorted(dic_keyword_name.keys(), key= lambda item: item.upper())

toc_all = '<a name="list"></a>'
toc_all += '<h3><b> Alphabetical list of all "How to ?" documentation topics.</b></h3>'

cur_let_all = 'A'
toc_all += "<p>"+cur_let_all+".&nbsp;\n"

# For each "topic" file
for keyword in list_keywords:
  topic_name=dic_keyword_name[keyword]
  topic_howto=dic_keyword_howto[keyword]
  while not (keyword.startswith(cur_let_all.lower()) or keyword.startswith(cur_let_all.upper())):
    cur_let_all = chr(ord(cur_let_all)+1)
    toc_all = toc_all + "<p>"+cur_let_all+".\n"
  toc_all = toc_all + "<br><a href=\"topic_"+ topic_name + ".html\">" + keyword + "</a> [How to "+topic_howto+"] &nbsp;&nbsp;\n"

#Generate a first version of the html file, in the order "header" ... up to the "end"
#Take the info from the component "default" if there is no information on the specific component provided in the yml file.
all_topics_html=""
for j in ["header","title","subtitle","copyright","links","toc_all","links","end"]:
  if j == "toc_all":
    all_topics_html += toc_all
  elif j == "subtitle":
    all_topics_html += "<h2>Complete list.</h2><hr>"
    all_topics_html += 'This document lists the names of all "How to ?" documentation topics for the abinit package.'
    all_topics_html += '<script type="text/javascript" src="../../js_files/generic_advice.js"> </script>'
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

activate_translation=0
if activate_translation==1:

  docs=yml_in["helps"]

  for i, doc_info in enumerate(docs):

    # Skip the default component
    name = doc_info.name
    if name=="default":
      break
    path_doc_html="users/"+name+"_help.html"
    path_doc_yml="users/origin_files/help_"+name+".yml"
    print("Read "+path_doc_html+" to build '"+path_doc_yml+"'... ",end="")

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
          # Correct the references to the other files in the topic directory (transient measure in case of the "what_ABINIT_does" file)
          string_old='href="generated_files/topic'
          string_new='href="../../topics/generated_files/topic'
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
            for i, var in enumerate(abinit_vars):
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
#            theory, files theorydoc_*
#            users,  files help_*
 
################################################################################

suppl_components={}
for list_infos in list_infos_dir:
  yml_files=list_infos["yml_files"]
  for yml_file in yml_files:
    if yml_file in ["lessons","theorydocs","helps"]:
      rc=assemble_html(yml_in[yml_file],suppl_components,list_infos["dir_name"],yml_file[:-1],allowed_link_seeds,backlinks)

################################################################################
################################################################################

# Come back to the bibliography

################################################################################
# Treat the links within the "introduction" of the acknowledgment component first.

backlink= ' &nbsp; <a href="../../bibliography/generated_files/acknowledgments.html">acknowledgments.html</a> &nbsp; ' 
for i, bibfile_info in enumerate(yml_in["bibfiles"]):
  if bibfile_info.name.strip()=="acknowledgments":
    bibfile_intro=bibfile_info.introduction
    bibfile_ack_intro = make_links(bibfile_intro,None,allowed_link_seeds,backlinks,backlink)

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
alphalinks+=" or <a href=#>top</a> \n \n"
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
  while ID[0]>cur_let:
    cur_let = chr(ord(cur_let)+1)
    bibliography_content+=('<a id="%s"></a>')%(cur_let)
    if cur_let==ID[0]:
      bibliography_content+=alphalinks+('<hr><hr><h2>%s</h2> \n \n')%(cur_let)
  bibliography_content+= ('<hr><a id="%s">[%s]</a> (<a href="./bibtex.html#%s">bibtex</a>)\n <br> %s \n') %(ID,ID,ID,reference_dic[ID])
  nlink=0
  for link in backlinksID:
    if len(link)!=0:
      if nlink==0: 
        bibliography_content+= "<br> Referred to in " 
        nlink=1
      bibliography_content+= link

# Open, write and close the txt file
file_txt = 'bibliography/generated_files/ordered_abiref.bib'
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

suppl={"introduction":bibfile_ack_intro}
suppl_components={"acknowledgments":suppl}
suppl={"content":bibtex_content}
suppl_components['bibtex']=suppl
suppl={"content":bibliography_content}
suppl_components['bibliography']=suppl

rc=assemble_html(yml_in["bibfiles"],suppl_components,"bibliography","",allowed_link_seeds,backlinks)

################################################################################

print("Work done !")

################################################################################
