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

def bibtex2html(str_input):
  """ Convert the bibtex notations to html notations inside the string str 
      The coding is often primitive and very specialized ... The goal is not to write a complete BibTex parser !
      If it does not work, modify the *.bibtex entry ...
  """

  str=str_input

  #Subscripts/superscripts without parentheses
  for i in string.digits:
    string_old='$_'+i+'$'
    string_new="<sub>"+i+"</sub>"
    str=str.replace(string_old,string_new)
    string_old='$^'+i+'$'
    string_new="<sup>"+i+"</sup>"
    str=str.replace(string_old,string_new)

  #Superscripts
  converted_str=convert_subscripts(str)
  str=convert_superscripts(converted_str)

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

  converted_str=convert_textit(str)
  str=suppress_parentheses(converted_str)
  str=str.replace("--","&ndash;")

  #Suppose remaining parentheses are present to avoid BibTex to switch automatically from uppercase to lowercase,
  #which will not happen in HTML...
  str=str.replace('"{','"')
  str=str.replace('}"','"')

  return str

################################################################################

def suppress_parentheses(text):
  def strip_text(mymatch):
    stripped_text = mymatch.group()[1:-1].strip()
    return stripped_text

  p=re.compile("\\{([a-zA-Z0-9_ */<>.|:+#@]*)\\}")
  if text is None:
    return ""
  stripped_text=p.sub(strip_text,text)

  return stripped_text

################################################################################

def convert_superscripts(text):
  def strip_text(mymatch):
    stripped_text = "<sup>"+mymatch.group()[2:-1].strip()+"</sup>"
    return stripped_text

  p=re.compile("\\$\\^\\{([a-zA-Z0-9_ */<>.|:+#@]*)\\}\\$")
  if text is None:
    return ""
  stripped_text=p.sub(strip_text,text)

  return stripped_text

################################################################################

def convert_subscripts(text):
  def strip_text(mymatch):
    stripped_text = "<sub>"+mymatch.group()[2:-1].strip()+"</sub>"
    return stripped_text

  p=re.compile("\\$_\\{([a-zA-Z0-9_ */<>.|:+#@]*)\\}\\$")
  if text is None:
    return ""
  stripped_text=p.sub(strip_text,text)

  return stripped_text


################################################################################

def convert_textit(text):
  def convert_text(mymatch):
    converted_text = "<i>"+mymatch.group()[8:-1].strip()+"</i>"
    return converted_text

  #Look for \textit{...}
  p=re.compile("\\\\textit\\{([a-zA-Z0-9_ */<>.|:+#@]*)\\}")
  if text is None:
    return ""
  converted_text=p.sub(convert_text,text)

  return converted_text

################################################################################
