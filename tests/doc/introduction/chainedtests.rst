=============
Chained Tests
=============

In some cases, one has to connect different runs 
A typical case is depicted in the figure below:

.. image:: dep.png

In this case, test t95.in produces a density file DEN that
is needed to start the second calculation in test t96.in.
The python code must be aware of this dependency so that 
t95.in and t96.in can be executed in the correct order. 

The TEST_INFO section defines the option test_chain 
that contains the name of the tests forming the chain.
The sequence is ordered in the sense that the python script will
execute the tests starting from the first item of the list.

.. note::
   All the input files belonging to a test chain must contain
   the test_chain option (obviously with the same list of tests)

During the execution of a test chain, one usually has to perform basic operations 
such as file renaming or file copying in order to connect the different steps 
(for example, we may want to rename the DEN file produced in t95.in so that 
the ouput density is automatically read in test t96.in).
The shell section contains two options that are used to specify the list of 
shell commands that are executed before and after running the test.
A typical example is shown below:

.. image:: test_chain.png


Note that we use a peculiar syntax for the shell commands

.. image:: rshell_syntax.png
