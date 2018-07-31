#!/usr/bin/env python
#
# Copyright (C) 2009-2018 ABINIT Group (Damien Caliste)
#
# This file is part of the ABINIT software package. For license information,
# please see the COPYING file in the top-level directory of the ABINIT source
# distribution.
#
# The purpose of this file is to check the documentation consistency for the
# input variables.
from __future__ import unicode_literals, division, print_function, absolute_import

import os
import re
try:
    import cPickle
except ImportError:
    import pickle as cPickle

test_dir       = ".."
top_dir        = os.path.join(test_dir,"..")
doc_dir        = os.path.join(top_dir,"doc")

# The documentation resources:
inp_dir        = os.path.join(doc_dir, "variables")
doc_index_file = os.path.join(inp_dir, "generated_files/varset_allvars.html")

# The source resource:
src_dir        = os.path.join(top_dir,"src")

# The dtset resource:
bindings_dir        = os.path.join(top_dir,"bindings")
parserdir = os.path.join(bindings_dir, "parser")

with open(os.path.join(parserdir, "dtset.pickle", "rb")) as fh:
    ab_dtset = cPickle.load(fh,"r")

# The documentation regexp.
re_index = re.compile(r"\<a\s+href\s*=\s*['\"](?P<file>var[a-z]+\.html)#(?P<var_link>[a-z0-9_]+)['\"]\s*\>\s*(?P<var_name>[a-z0-9_,\s]+)\s*\</a\>", re.IGNORECASE + re.MULTILINE)
re_title = re.compile(r"\<font\s+id\s*=\s*['\"]title['\"]\>\s*\<a\s+name\s*=\s*['\"](?P<var_id>[a-z0-9_]+)['\"]\>\s*(?P<var_name>[a-z0-9_,\s]+)\s*\</a\s*\>\s*\</font\s*\>", re.IGNORECASE + re.MULTILINE + re.DOTALL)
# The intagm routine regexp.
re_intagm = re.compile(r"^ *call +intagm *\(([a-zA-Z, 0-9():'%_*+-/]+) *&?\n? *&? *([a-zA-Z, 0-9():'%_*+-/]*)\) *$", re.MULTILINE + re.DOTALL + re.IGNORECASE)

def src_parse_file(arg, dirname, fnames):
  (re_intagm, ab_vars) = arg
  for f in fnames:
    filename = os.path.join(dirname, f)
    if (os.path.isfile(filename) and filename.endswith(".F90")):
      with open(os.path.join(dirname, f), "rt") as tmpFile:
        src_file = tmpFile.read()

      code = re_intagm.split(src_file)
      for i in range((len(code) - 1) / 3):
        if (re.match('^ *[^!] *call +intagm', code[i * 3]) is not None):
          print(code[3 * i])
          print(filename)
          raise ValueError
        (token, tread) = (code[i * 3 + 1] + code[i * 3 + 2]).split(",")[6:8]
        re_set = re.compile(re.escape(token) + " *= *['\"]([a-z0-9][a-z0-9_]*)\s*['\"]", re.IGNORECASE + re.MULTILINE)
        var = None
        for j in range(i, -1, -1):
          res = re_set.search(code[j * 3])
          if (res is not None):
            var = res.group(1)
            break
        if (var is None):
          print(filename)
          print(token)
          raise ValueError
        # Test if we don't leave_new after read of this variable.
        if (re.match(r"^\s*if\s*\(\s*" + re.escape(tread) + \
                     r"\s*==\s*1\s*\)\s*then\s*write.+" + \
                     r"call\s+wrtout.+" + \
                     r"call\s+leave_new.+" + \
                     r"end\s*if", code[3 * i + 3], \
                     re.IGNORECASE + re.MULTILINE + re.DOTALL) is None):
          ab_vars[var] = [filename]

def doc_parse_index(filename, defs, ab_stat_links):
  # Load the documentation into memory for later search.
  with open(os.path.join(inp_dir, filename), "rt") as tmpFile:
    doc_index = tmpFile.read()

  # Put all the matching elements into an iterator
  keysIter = re_index.finditer(doc_index)

  # Consistency check for the links.
  links = {}
  for var in keysIter:
    ele = var.groupdict()
    names = ele['var_name'].split(',')
    for name in names:
      name = name.strip()
      links[name] = [False]
      ab_stat_links[0] += 1
      if name in defs and \
           (not(defs[name][2] == ele['var_link']) or not(defs[name][1] == ele['file'])):
        print("Warning: variable '%11s' as a wrong link '%s#%s' in the HTML file '%s' (should be '%s#%s')." % \
            (name, ele['file'], ele['var_link'], filename, defs[name][1], defs[name][2]))
        ab_stat_links[1] += 1
  return links

def doc_parse_file(arg, dirname, fnames):
  (re_title, ab_defs) = arg
  for f in fnames:
    filename = os.path.join(dirname, f)
    if (os.path.isfile(filename) and f.startswith("var")):
      with open(os.path.join(dirname, f), "r") as tmpFile:
        doc_src = tmpFile.read()

      # Put all the matching elements into an iterator
      keysIter = re_title.finditer(doc_src)

      # Store all the definitions.
      for var in keysIter:
        ele = var.groupdict()
        names = ele['var_name'].split(',')
        for name in names:
          name = name.strip()
          if name in ab_defs:
            print("Warning: variable '%11s' from '%s' has already a definition in '%s'." % \
                  (name, f, ab_defs[name][1]))
          ab_defs[name] = [False, f, ele['var_id']]
        if (len(names) == 1) and not(ele['var_id'] == ele['var_name']):
          print("Warning: variable '%11s' as a wrong name anchor '%s' in the HTML file '%s'." % \
                (ele['var_name'], ele['var_id'], f))
          ab_defs[name][2] = name

ab_stat_links = [0, 0]

print("Dtset structure:")
print("===============")
print("There are %d variables in the dataset structure.\n" % len(ab_dtset))

print("var*.html files:")
print("===============")
# Parse all var*.html files.
ab_defs = {}
os.path.walk(inp_dir, doc_parse_file, (re_title, ab_defs))
# Check that links to vars are consistent in all var* files.
for f in os.listdir(inp_dir):
  if (f.startswith("var")):
    doc_parse_index(os.path.basename(f), ab_defs, ab_stat_links)
print("There are %d variables defined in the var*.html files.\n" % len(ab_defs))

print("HTML index file:")
print("===============")
keys = doc_parse_index("generated_files/varset_allvars.html", ab_defs, ab_stat_links)
print("There are %d variables in the HTML index file.\n" % len(keys))

print("Intagm calls:")
print("============")
ab_vars = {}
os.path.walk(src_dir, src_parse_file, (re_intagm, ab_vars))
print("There are %d variables read from intagm().\n" % len(ab_vars))
#print("ab_vars.keys", ab_vars.keys())

# Check that all ab variables have been defined in the index file.
ab_no_def = 0
for (ab_var, ab_att) in ab_vars.items():
  no_index = True
  no_def = True
  if (ab_var in keys):
    keys[ab_var][0] = True
    no_index = False
  if (ab_var in ab_defs):
    ab_defs[ab_var][0] = True
    no_def = False
  if no_index and not(no_def):
    print("Warning: ABI var. '%11s' (%s) is missing from the HTML index file."% (ab_var, ab_att[0]))
    ab_no_def += 1
  elif not(no_index) and no_def:
    print("Warning: ABI var. '%11s' (%s) is not defined in any var*.html file." % (ab_var, ab_att[0]))
    ab_no_def += 1
  elif no_index and no_def:
    print("Warning: ABI var. '%11s' (%s) is not defined anywhere."% (ab_var, ab_att[0]))
    ab_no_def += 1
# Check that all dtset variables have been defined in the index file.
for dt_var in list(ab_dtset.keys()):
  if (dt_var in keys):
    keys[dt_var][0] = True
  if (dt_var in ab_defs):
    ab_defs[dt_var][0] = True
  
# Check that there is no unused variable in the index file.
for (key, var) in keys.items():
  if (var[0] is False):
    print("Warning: variable '%11s' is declared in the index file but not used." % key)
for (key, var) in ab_defs.items():
  if (var[0] is False):
    print("Warning: variable '%11s' (%s) is declared but not used in the code." % (key, var[1]))

print()
print("Summary:")
print("=======")
print(" - %d/%d (%4.2f%%) of documented input variables." % ((len(ab_vars) - ab_no_def), len(ab_vars), (float(len(ab_vars) - ab_no_def) / float(len(ab_vars)) * 100.)))
print(" - %d/%d (%4.2f%%) of broken links." % (ab_stat_links[1], ab_stat_links[0], (float(ab_stat_links[1]) / float(ab_stat_links[0]) * 100.)))
