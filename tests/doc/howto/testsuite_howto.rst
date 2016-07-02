====================
How to run the tests
====================

The script abinit/tests/runtests.py provides a user-friendly interface
that allows the user to select tests by keywords, authors, suite name, input variables etc. 
It also provides options for controlling the number of MPI processes, 
the number of OpenMP threads and the number of python threads (task parallelism).

The script syntax is::

    runtests.py [suite_args] [options]. 

where suite_args is a list of suite names that can be optionally selected (sliced) using python syntax

For example::

    runtests v3[:4] v4[45:] v5[3] 

execute all the tests in v3 from 0 up to 4 (excluded), all the tests in v4 starting from 45 and v5/t3

.. warning::
    If suite_args is not given, the script will execute the entire list of tests.

The most useful options are:

- -j:
    specifies the number of python threads (task parallelism)

- -n:
    specifies the number of MPI nodes

- -o:
    specifies the number of OMP threads

- -k:
    selects tests by keywords

- -h:
    for help + list of examples

.. tip::
    For the complete list of options supported, use runtests -h

=========================
How to run tests with MPI
=========================

In the simplest case, one can run MPI tests by just using::

    runtests.py -n 2 paral

where n specifies the number of MPI processes.
runtests.py employs very simple rules for guessing the MPI implementation that must be invoked. 
In particular, it assumes that:

    1. The environment is properly setup (PATH, LD_LIBRARY_PATH)
    
    2. MPI Demons (e.g. mpd) are already running in the background 
       (you have to initialize the demon manually before launching the tests)

    3. It uses the mpirunner in PATH (default: mpirun) and assumes 
       the standard syntax for specifying the number of processors (i.e. mpirun -n 4 exe < stdin > stdout)

If you need to override the default setting, you can pass options to the script
via a configuration file in the INI format (-c options)::

    runtests.py -n 2 -c mpi.cfg

    $ cat mpi.cfg

        [mpi]
          mpi_prefix = /usr/local/openmpi-gcc47/
          mpirun_np = %(mpi_prefix)s/bin/mpirun -np

mpirun_np is the string with the path to the mpirunner, followed by any additional option you may want
to use. The option used to specifing the number of MPI nodes must be added at the end of the string.

.. tip::
    For more examples, consult the configuration files in abinit/tests/mpi_cfg.

=========================================
How to run testbot.py on a buildbot slave 
=========================================

This section is intended for developers who might want to run interactively 
the entire test suite on a slave builder (example: one or more tests 
are failing on a slave. You have already modified the source code to fix the problem
and you want avoid to launch a full builbot build just to check that your changes are OK)

In this case, one can re-rerun the entire test suite (or part of it) by just executing the
following two steps:

    $ cd abinit/tests
    $ testbot.py

The script testbot.py reads the configuration file testbot.cfg (already present in the working directory), 
runs the entires set of tests and produces the final report.

Note that one modify the configuration options defined in testbot.cfg in order to speed-up the execution of the tests.
In particular one can use the options:

    # with_tdirs = list of directories to execute
    # without_tdirs = list of directories that will be excluded

=====================
How to add a new test
=====================

Let's assume you want to add a new test tfoo_1.in to the directory v5.
This test will be hence located in the subsuite foo within the suite v5 
so that we can easily run it with runtests by using::
    
    runtests foo[1] 

In order to add the new test, one has to follow the following steps:

    1. Add the name of the input file to the python list inp_files 
       defined in the file v5/__init__.py

    2. Register the name of the subsuite in the list subsuites 
       defined v5/__init__.py (this step can be skipped if the 
       test does not belong to a subsuite)

    3. Add the input files to v5/Input and the reference files to v5/Refs

    4. Run teslint.py to make sure that your modifications do not break
       the basic rules and the conventions assumed by the python code.
   
.. warning::

       The code uses a Pickle database (test_suite.cpkl) to store the set of objects
       representing the tests of the test suite.
       Remember to regenerate the database after any change to the TEST_INFO section
       or any modification of the configuration parameters of the test suite (__init__.py files).
       runtests.py provides the handy option -r (regenerate) to automatically regenerate the database
       before running the tests.

===========================
How to add a chain of tests
===========================
TODO

==========================
How to add a parallel test
==========================
TODO

=======================================
How to add support for a new executable
=======================================

Let's assume that you have a new executable named foo.x
and you want to change the python code so that the results of foo.x are automatically tested.
What are the modifications needed to integrate foo.x in the ABINIT test suite?

First you have to provide an input file for foo.x with a
<TEST_INFO> section that provides ALL the information needed to run the executable. 
If foo.x requires some kind of input data or extra rules that are not supported, 
you will have to modify the parser so that the new options are stored in BaseTest

Then you have to tell the code how to construct the standard input that will be passed to foo.x. 
This is the most complicated part as it requires some understanding of the internal implementation.

The code uses three different objects:
 
    1. BaseTest
    2. ChainOfTests
    3. TestSuite

to represent the tests of the testsuite.

BaseTest is the base class that provides methods to run and analyze the results
of the test (you should try to reuse this piece of code as much as possible).

ChainOfTests is a list of BaseTest instances that cannot be executed separately. 
Typical example: t1.in produces an output result that is used as the input of tests t2.in.

TestSuite is a list of (BaseTest, ChainOfTests) instances and provides
user-friendly methods to extract tests according to some rule (e.g. keywords) 
and to run these tests with Python threads.

The BaseTest class provides the method make_stdin that returns 
a string containing the standard input that should be passed to the Fortran executable.
Since each executable has its own format for the input file, the BaseTest is not able to handle all the different cases.
In many cases, you only have to replace the make_sdtin method of BaseTest 
so that the appropriate standard input is constructed from the values specified in the <TEST_INFO> section.  

