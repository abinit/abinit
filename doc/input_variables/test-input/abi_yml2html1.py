#! /usr/bin/env python

import yaml
import re
import argparse
from variables import *
debug = 0

# Path relative from HTML files
js_path = ""
users_path = "../../../users/"
section_path="../../input_variables/html_automatically_generated/"

file='abinit_vars.yml'

parser = argparse.ArgumentParser(description='Tool for eigenvalue analysis')
parser.add_argument('-f','--file',help='YML file to be read')
args = parser.parse_args()
args_dict = vars(args)
if args_dict['file']:
  file = args_dict['file']

print("Will use "+str(file)+" as input file")

with open(file, 'r') as f:
    variables = yaml.load(f);

##output = ''

with open('temp_allvariables.html') as f:
    header_all = f.read()

with open('temp_specs.html') as f:
    header_specs = f.read()

header_all = header_all.replace("__JS_PATH__",js_path)

all_contents = dict()
all_vars = dict()

list_all_vars = dict()

for var in variables:
  list_all_vars[var.varname] = var.section

# Recompute groups
for i, var in enumerate(variables):
 if debug==1 :
   print(var)
 varname = var.varname
 if var.characteristics is not None and '[[INTERNAL_ONLY]]' in var.characteristics:
   varname = '%'+varname
 section = var.section
 print(all_contents)
 if section not in all_contents.keys():
   all_contents[section] = "<br><br><br><br><hr>\n"
   all_vars[section] = []

 all_vars[section].append([var.varname,var.definition])

for section, content in all_contents.items():
 file_cur = ''+section+'.html'
 f_cur = open(file_cur,'w')

 with open('temp_'+section+'.html') as f:
    header_varX = f.read()

 cur_header_varX = header_varX.replace("__JS_PATH__",js_path)
 f_cur.write(cur_header_varX)

 for varname,defi in all_vars[section]:
  curlink = " <a href=\""+section_path+section+".html#"+varname+"\">"+varname+"</a>&nbsp;&nbsp;\n"
  f_cur.write(curlink)

 f_cur.write("\n")
 content += "<script type=\"text/javascript\" src=\""+js_path+"list_internal_links.js\"> </script>\n\n"
 content += "</body>\n"
 content += "</html>"
 f_cur.write(content)
 f_cur.write("\n")
 f_cur.close()


