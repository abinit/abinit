#!/usr/bin/env python
from __future__ import print_function, division, absolute_import #, unicode_literals

import os
import sys
import imp

from os.path import join as pj, abspath as absp, exists as pexists, basename, splitext, isfile

#try:
#  import tests as abitests
#except ImportError:
# Add the directory [...]/abinit/tests to $PYTHONPATH
pack_dir, x = os.path.split(absp(__file__))
pack_dir, x = os.path.split(pack_dir)
sys.path.insert(0,pack_dir)
import tests as abitests
del pack_dir, x

def tests_from_pymod(pymod_path, abenv):

  # Import the module
  mod_name = basename(pymod_path).split(".py")[0] 
  module = imp.load_source(mod_name, pymod_path)

  # Inspect the module and build the list of tests.
  tests = list()

  test_gen = getattr(module,"abinit_test_generator", None)
  if test_gen is not None:
    tests.append( PythonTest(abenv, test_gen()) )

  suite_gen = getattr(module, "abinit_suite_generator", None)
  if suite_gen is not None:
    for test_gen in suite_gen(): 
       tests.append( PythonTest(abenv, test_gen) )

  assert tests
  return tests

class PythonTest(object):
  _attrbs = [
    "test_func",
    #"pyscript"
    #"name",
  ]

  def __init__(self, abenv, dictionary):
    self.abenv = abenv

    for k in PythonTest._attrbs:
       self.__dict__[k] = dictionary.get(k, None)

    assert self.test_func is not None
    self.info = self.test_func.__doc__
    assert self.info is not None

    # TODO: Import conditions.
      
  def __str__(self):
    return "\n".join( [ str(k) + " : " + str(v) for k,v in self.__dict__.items()] )

  def run(self, workdir=None):
    if workdir is None: workdir = os.path.abspath(os.curdir) 
    self.workdir = workdir

    print(self.info +"...",)
    #stderr_fname = pj(workdir, self.name + ".stderr")
    #stdout_fname = pj(workdir, self.name + ".stdout")

    #self.stderr = open(stderr_fname, "w")
    #self.stdout = open(stdout_fname, "w")
    from StringIO import StringIO 
    self.stdout, self.stderr = StringIO(), StringIO()

    sys.stdout, sys.stderr = self.stdout, self.stderr

    try:
      exit_status = self.test_func(self.abenv)
    except:
      import traceback
      var = traceback.format_exc()
      #print(self.name + ": " + var)
      print(var)
      exit_status = 99

    self.stdout.seek(0)
    self.stderr.seek(0)

    #self.stderr.close()
    #self.stdin.close()

    # Reinstate stdout and stderr.
    sys.stdout, sys.stderr = sys.__stdout__ , sys.__stderr__ 

    #print("stdout:")
    #print(self.stdout.read())
    #print( "stderr:")
    #print( self.stderr.read())
    #print( "done")

    if exit_status != 0: # Check reference files.
      print(" (Comparing with reference files)")
      #FIXME: I need to change the name of the reference files.
      #print self.stdout.read()
      #print self.stderr.read()

      ref_dir = self.abenv.apath_of("tests/abirules/Refs")
      #for idx, ext in enumerate(["", ".stdout", ".stderr"]):
      #for idx, ext in enumerate([""]):
      for idx, ext in enumerate([]):
        #ref_file = pj(ref_dir, self.name + ext)
        ref_file = pj(ref_dir, "t01.out")
        if not isfile(ref_file): continue
        fh = open(ref_file, "r")
        ref_lines = fh.readlines()
        fh.close()
        if idx == 0: from_lines = self.stdout.readlines()
        if idx == 1: from_lines = self.stderr.readlines()
        import difflib
        diff = difflib.unified_diff(from_lines, ref_lines) #, fromfile, tofile, fromdate, todate, n=n)
        diff = list(diff)
        for l in diff: print(l)
        exit_status = 0
        if len(diff) != 0: exit_status = 1
        
    self.exit_status = exit_status

    if exit_status != 0: 
      print("FAILED")
    else:
      print("SUCCESS")

    return self.exit_status

def run_abirules_suite():
  abenv = abitests.abenv
  import abirules 
  script_dir = abenv.apath_of("abichecks/scripts/")
  scripts = [ absp(pj(script_dir,s)) for s in abirules.pyscripts ]

  exit_status = 0
  for script in scripts:
    tests = tests_from_pymod(script, abenv)
    for test in tests:
      retval = test.run()
      if retval != 0 : exit_status = retval
  return exit_status

def run_buildsys_suite():
  abenv = abitests.abenv

  import buildsys
  script_dir = abenv.apath_of("abichecks/scripts/")
  scripts = [ absp(pj(script_dir,s)) for s in buildsys.pyscripts ]

  exit_status = 0
  for script in scripts:
    tests = tests_from_pymod(script, abenv)
    for test in tests:
      retval = test.run()
      if retval != 0 : exit_status = retval
  return exit_status

if __name__ == "__main__":
 
  #print 30*"=" + " buildsys_tests " + 30*"="
  #run_buildsys_suite()

  print(30*"=" + " abirules_tests " + 30*"=")
  run_abirules_suite()
