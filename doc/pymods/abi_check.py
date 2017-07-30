#!/usr/bin/env python
from __future__ import print_function

import sys

from variables import *

try:
  import yaml
except ImportError:
  raise ImportError("yaml package is not available. Install it with `pip install pyyaml`")


def main():
  path = '../input_variables/origin_files/'
  path_topics = '../topics/origin_files/'
  path_abinit_vars = path+'abinit_vars.yml'
  with open(path+'characteristics.yml', 'r') as f:
    characteristics = yaml.load(f)
  with open(path+'varsets.yml', 'r') as f:
    list_varsets = yaml.load(f)
    varset_names=[]
    for varset in list_varsets:
      varset_names.append(varset.name)
  with open(path_topics+'list_of_topics.yml', 'r') as f:
    topics = yaml.load(f)
    for i,item in enumerate(topics) :
       topics[i]=item.strip()
  with open(path_topics+'list_tribes.yml', 'r') as f:
    tribes = yaml.load(f)
    tribe_names=[]
    for item in tribes:
      tribe_names.append(item[0].strip())

  retcode = 0

  #Detect incorrect YAML syntax (at least, the ones not compatible with the yaml python module ...) by brute force
  #specific means
  with open(path_abinit_vars,'r') as f:
    abinit_vars=f.readlines()
    for i,line in enumerate(abinit_vars):
      if ("abivarname"==line.lstrip()[:10] or
          "characteristics"==line.lstrip()[:15] or
          "commentdefault"==line.lstrip()[:14] or
          "commentdims"==line.lstrip()[:11] or
          "defaultval"==line.lstrip()[:10] or
          "dimensions"==line.lstrip()[:10] or
          "excludes"==line.lstrip()[:8] or
          "executables"==line.lstrip()[:11] or
          "mnemonics"==line.lstrip()[:10] or
          "requires"==line.lstrip()[:9] or
          "text"==line.lstrip()[:4] or
          "topics"==line.lstrip()[:6] or
          "varset"==line.lstrip()[:6] or
          "vartype"==line.lstrip()[:7]):
        line_split=line.strip().split(":",1)
        if len(line_split)>1:

          if line_split[1] !=""  and line_split[1][0] != " ":
            print('FAIL: detected a non-blank character after the first colon in the line number %s :'%(i+1))
            print(line)
            print('This is not allowed by YAML syntax.')
            retcode+=1

          line1=line_split[1].strip()
          if len(line1)>0:
            line1_split=line1.strip().split(":",1)
            if len(line1_split)>1:
              if not  ( "'" in line1_split[0] or
                        '"' in line1_split[0] or  
                        "(" in line1_split[0] or 
                        "[" in line1_split[0] ):
                print('FAIL: detected a second colon in the line %s :'%(i+1))
                print(line)
                print(line1)
                print("line1_split:",line1_split) 
                print('This is not allowed by YAML syntax.')
                retcode+=1

  if retcode>0:
    sys.exit()

  #Load abinit_vars.yml using YAML syntax.
  with open(path_abinit_vars, 'r') as f:
    variables = yaml.load(f)

  for var in variables:

    if var.abivarname is None:
      print('Variable %s has no name. This is forbidden' % str(var))
      retcode += 1
      continue

    abivarname = var.abivarname

    if var.vartype is None:
      print('FAIL: ', abivarname, ' does not have a vartype')
      retcode += 1
    elif not var.vartype in ["integer", "real", "string"]:
      print('FAIL: ',abivarname,' has vartype ',var.vartype,' not in the list ["integer", "real", "string"].')
      retcode += 1

    if var.topics is None:
      print('FAIL: ', abivarname, ' does not have at least one topic and the associated tribe')
      retcode += 1
    else:
      topics_name_tribe = var.topics.split(',')
      for topic_name_tribe in topics_name_tribe:
        name_tribe = topic_name_tribe.split('_')
        if not name_tribe[0].strip() in topics:
          print('FAIL: ', abivarname, ' delivers topicname_tribe ',name_tribe,
                ' with topicname ',name_tribe[0].strip(),' that does not belong to the allowed list')
          retcode += 1
        if not name_tribe[1].strip() in tribe_names:
          print('FAIL: ', abivarname, ' delivers topicname_tribe ',name_tribe,
                ' with tribe ',name_tribe[1].strip(),' that does not belong to the allowed list')
          retcode += 1

    if var.characteristics is not None:
      if not isinstance(var.characteristics, list):
        print('FAIL: the field characteristics of ', abivarname, ' is not a list')
        retcode += 1
      else:
        for cat in var.characteristics:
          if cat.replace("[[", "").replace("]]", "") not in characteristics:
            print('FAIL: the characteristics ', cat, ' of ', abivarname, ' is not valid')
            retcode += 1

    if var.dimensions is not None:
      if var.dimensions != "scalar":
        if not isinstance(var.dimensions, list) \
              and not isinstance(var.dimensions, ValueWithConditions):
          print('FAIL: the field dimensions of ', abivarname, ' is not a list neither a valuewithconditions')
          retcode += 1
    else:
      print('FAIL: ', abivarname, ' does not have a dimension. If it is a "scalar", it must be declared so.')
      retcode += 1

    if var.varset is None:
      print('FAIL: ', abivarname, ' does not have a varset')
      retcode += 1
    else:
      if not isinstance(var.varset, str) or var.varset not in varset_names:
        print('FAIL: the field varset of ', abivarname, ' should be one of the valid varsets')
        retcode += 1

    varname_split=abivarname.split("@",1)
    if len(varname_split[0])>20:
      print('WARNING: len of principal name of ',abivarname,', namely, ',len(varname_split[0]),', is longer than 20 characters.')

  if retcode != 0:
    print('Found ',retcode,' FAIL.')
    print('Number of variables in ABINIT :',len(variables))

  return retcode


if __name__ == "__main__":
  sys.exit(main())
