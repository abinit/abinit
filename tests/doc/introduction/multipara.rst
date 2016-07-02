====================
Multi-parallel Tests
====================

A multi-parallel test is a test that can be executed in parallel
with different numbers of MPI processes (e.g. the tests in tests/paral).

The ABINIT test suite is designed so that it is possible to define a single input file
from which multiple parallel tests are automatically generated.
In what follows, we describe how to specify a multiparallel test in a TEST_INFO section. 

The most important option is nprocs_to_tests (defined in the [mpi] section) 
that lists the number of MPI processes to use for the different tests. 
max_nprocs is the maximum number of MPI processes 
that are supported by the calculation (just set it to the maximum value present in nprocs_to_test)

Then for each possible number of MPI processes we have a section whose name 
is constructed with the rule: NCPU_#nproc

A simplified example of TEST_INFO section for a multi-parallel test is reported below

.. image:: multipara.png


==========================================================
A more complicate example: a chain of multi parallel tests
==========================================================

.. warning::
   You cannot change the number of MPI processes inside a test chain
