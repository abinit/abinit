======
Design
======

The ABINIT test suite is a collection of input and output files organized in a hierarchical structure of directories. 
Each directory contains tests that share some kind of property such as 
the ABINIT major version in which the test was added or the external libraries whose presence 
is mandatory for the execution of the test.

The purpose of the test suite is two-fold: developers use the test farm
to validate their code before pushing their developments to public whereas 
users use the test suite to validate the executables used for their production calculations.

The code that drives the execution of the tests and the validation of the final results 
in mainly written in python (version >= 2.4 is required, not compatible with 3.0) 
and relies on the following principles:

    #.  Each test is self-describing, in the sense that the input file 
        provides all the information needed to run the test and to analyze the final results.

    #.  Each test has some kind of keyword so that it is possible to 
        select and run only the tests that contain these tags. 

    #.  It should be possible to run easily the entire suite with MPI and an arbitrary 
        number of processors (not only the tests in paral)

    #.  The execution of the tests can be distributed among python threads in order 
        to accelerate the execution, in particular on multi-core architectures.

===================
Configuration files
===================

The configuration parameters of the test suite are defined in the __init__.py
files located in the test suite directories.
These files are standard python files, the parameters are expressed with the standard python syntax.

tests/__init__.py is the top level configuration file that defines the list 
of directories containing ABINIT tests::

    cat tests/__init__.py

    testsuite_dirs = [
      "atompaw",
      "bigdft",
      "built-in",
      "etsf_io",
      "fast",
      #"cpu",  <<< disabled
      ...
    ]

To add a new directory to the test suite, it is sufficient to add its name to the python list.

Each directory listed in testsuite_dirs contains an additional __init__.py file 
that specifies the list of input files associated to the tests and (optional) global attributes::

    cat tests/libxc/__init__.py

    # CPP variables 
    need_cpp_vars = [
      "HAVE_LIBXC",
    ]

    # List of input files
    inp_files = [
    "t00.in",
    "t01.in",
    "t02.in",
    "t03.in",
    "-t11.in", # Disabled
     ...
    ]

inp_files is the list with the names of the input files contained in the Input directory.
One can prepend a dash to the name of file (e.g. “-t11.in”) to signal to the python code 
that this test has been disabled and it should not be executed on the test farm.

The variable need_cpp_vars is a list with the CPP variables that must be defined 
in the include file config.h in order to enable the tests in this directory.
In this case, for example, the libxc tests are executed only if HAVE_LIBXC is defined in config.h
Note that one can prepend the character ! to the name of the variable to specify that the tests 
should not be executed if the variable is not defined in the build.

Each input file contains all the information needed to run the test and to analyze the 
final results. These meta-variables are gathered in the so-called TEST_INFO section
that is discussed in more detail in ...

===========
Conventions
===========

This paragraph summarizes the conventions used for the names of the tests.
The ABINIT tests are grouped in suites whose name is given by the name of the directory
(as specified in tests/__init__.py). Each suite can be optionally divided in sub-suites
whose names must be registered in dirname/__init__.py

.. note::
   The name of the suite/subsuite must be unique. 

To each test is therefore assigned a suite, a sub-suite and an integer number (>=0) 

The organization in terms of suites/subsuites allows one to select easily the tests
with the command line interface runtests.py discussed in ...

The name of the input file must match one of the two possible regular expressions:

    1. t[integer].in  e.g. t11.in

    2. t[subsuitename]_[integer].in  e.g. tgw1_1.in

Pattern 2) is used, for example, inside the tutorial directory to 
group the tests associated to a particular lesson.
The name of sub-suite can be passed as argument to runtests.py to specify
that only the tests in that particular subsuite should be executed.

For example, the command::
  
   runtests.py gw1[1:3]

runs the tests tgw1_1.in, tgw1_2.in.

.. warning::
    Each test has a standard input, a standard output and a standard error 
    whose names are tgw1_1.stdin, tgw1.stdout, and tgw1.stderr respectively.
    These names are constructed automatically by the python code hence 
    developers should avoid creating similar file names in the Fortran executables at runtime.
