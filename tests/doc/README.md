abinit test suite documentation
===============================

This is the top level build directory for the abinit testsuite documentation.  
All of the documentation is written using sphinx, a python documentation system 
built on top of ReST. This directory contains

* index.rst - the top level include document 

* introduction - Documents describing the design of the test suite

* howto - documentation with HOWTOs

* api - placeholders to automatically generate the api documentation

* conf.py - the sphinx configuration

* _static - used by the sphinx build system

* _templates - used by the sphinx build system

To build the HTML documentation, install sphinx (1.0 or greater
required), then type "make -f _Makefile html" in this directory.  
The documentation is produced in _build (see _build/index.html).
