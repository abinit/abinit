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
  with open(path+'abinit_vars.yml', 'r') as f:
    variables = yaml.load(f)
  with open(path+'characteristics.yml', 'r') as f:
    characteristics = yaml.load(f)
  with open(path+'varfiles.yml', 'r') as f:
    list_varfiles = yaml.load(f)
    varfile_names=[]
    for varfile in list_varfiles:
      varfile_names.append(varfile.name)
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
  for var in variables:

    if var.abivarname is None:
      print('Variable %s has no name. This is forbidden' % str(var))
      retcode += 1
      continue

    abivarname = var.abivarname

    if var.vartype is None:
      print('FAIL: ', abivarname, ' does not have a vartype')
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

    if var.varfile is None:
      print('FAIL: ', abivarname, ' does not have a varfile')
      retcode += 1
    else:
      if not isinstance(var.varfile, str) or var.varfile not in varfile_names:
        print('FAIL: the field varfile of ', abivarname, ' should be one of the valid varfiles')
        retcode += 1
 
  if retcode != 0:
    print('Found ',retcode,' FAIL.')
    print('Number of variables in ABINIT :',len(variables))

  return retcode


if __name__ == "__main__":
  sys.exit(main())
