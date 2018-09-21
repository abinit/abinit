---
authors: MG, XG
---

# HowTo guide for developers

This page is intended as a quick reference to solve problems commonly encountered
when developing new features in Abinit.


<!--
## Checklist

Checklist

1. Introduction.
2. Check list

*******************************************************************

1. Introduction

A few days before the developers sent their contributions
for the merge, it is worth that they examine whether they
have not forgotten to do "something". 
This is the aim of this checklist.

2. Check list

Have you mentioned your initials (or full name) in the
copyright of the routine you have modified ?

Have you set up one (or a few) test cases, so as to provide
the input and output files as references, to check that
the merge of your routines has been succesful?
The CPU time of this test should exceed 30secs on a PC at 400MHz,
or 15 secs on a better workstation.
-->

## How to add a new input variable 

!!! tip

    I assume you already know how to use [ctags](http://ctags.sourceforge.net/)
    to find routines, functions and datatypes in the source tree
    (you don't use grep to browse the code, do you?)
    so I'm not going to provide the full path to the F90 files as one can easily open the
    declaration of `dataset_type` with:

        cd ~abinit/src && ctags -R && vi -t dataset_type

    See [here](https://andrew.stwrt.ca/posts/vim-ctags/) for more tips.


Let's focus on the procedure required to add a new Abinit variable.
To make things as simple as possible, we neglect the case of dimensions such as *nkpt* or *nsym*
whose value may depend on the dataset.

To add a new variables follow the below steps:

- Add the new variable to `dataset_type`. Remember that the name cannot end with a digit (multidataset syntax)

- The default value can be specified in two different ways:

    * in the **declaration** of the Fortran type if the size is known at compile time 
      and the initial value does not depend on other variables.
    * in the **indefo** routine if the value must be computed at runtime.

- Add the name of the new variable to **chkvars**.

- Add a new section to **dtset_copy** to copy the new variable (use **alloc_copy** if allocatable).

- If you need an **allocatable entity**, remember to **deallocate** memory in **dtset_free**.

- Read the variable in the **invars2** (if it is not a basic dimension).

- Change one of the outvars routines (**outvar_a_h**, **outvar_i_n**, **outvar_o_z**) to print the variable 
  according to the first letter of the new variable

- The logic for checking the consistency of input variables goes to **chkinp**.
  Use the routines *chkint_eq*, *chkint_ne*, *chkint_ge*, *chkint_le*, *chkdpr*.

Finally,

    make clean && make -j8

since you *broke* the [ABI](https://en.wikipedia.org/wiki/Application_binary_interface) 
of a public datastructure and all the object files depending on this datastructure must be recompiled 
(if you are developing a library, you should release a new major version!)

No, it's not a typo, ABIs and APIs are different concepts!
From this detailed answer on [stackoverflow](https://stackoverflow.com/questions/2171177/what-is-an-application-binary-interface-abi)

>   If you expand, say, a 16-bit data structure field into a 32-bit field, then already-compiled code 
    that uses that data structure will not be accessing that field (or any following it) correctly. 
    Accessing data structure members gets converted into memory addresses and offsets during compilation 
    and if the data structure changes, then these offsets will not point to what the code is expecting 
    them to point to and the results are unpredictable at best.


For the treatment of dimensions see **invars0**, **invars1m**

## How to add a new test in the test suite?

The following information complements the [testsuite documentation](testsuite_howto).

In order to introduce a test, one needs to:

-  Provide a new input file in e.g. **tests/v8/Input**
   (or modify an existing one, possibly in another directory, but please do not suppress the existing capability testing!)

-  Provide a new reference file in the corresponding **tests/v8/Refs**

-  Insert the test in the list of tests to be done, by adding its name in **tests/v8/\_\_init\_\_.py**.

-  Document the test by adding a commented section inside the input file 
   (edit an existing input file to follow its style)

- Declare the pseudopotentials, the flow of actions, the files to be analyzed, the tolerances for the test, 
  inside this documentation. For the tolerances, start with 
  
        tolnlines = 0, tolabs = 0.000e+00, tolrel = 0.000e+00
    
  The different fields control how strictly the test results will be analysed:

  *  the maximum floating point difference found must be less than tolabs,
  *  the maximum relative difference less than tolrel, for each line individually,
  *  **tolnlines** is the maximal number of lines found to be different between the output and reference file 
     (within another tolerance that can be tuned by the option **opt=-medium**, **opt=-easy** or **opt=-ridiculous**,
     see the added section in other input files).

The scripts *tests/Scripts/fldiff.pl* and *tests/Scripts/reportdiff.pl* analyze the output.

If this procedure fails, contact the code maintainers in order adjust the test case. 
This might mean modifying the tolerances files. 
Unless you are really expert, let the maintainer do the final adjustment.

!!! important

    Please, try to keep preferably the total time per test case to less than 10 seconds.  
    30 seconds is a maximum for most test cases. Going being this needs exceptional reasons.

Supposing now that you have introduced a new test case. 
It can be used, with the other tests of the ABINIT test suite, through different channels:

*  these tests can be triggered on the test farm Web Page.

*  locally (on your machine), the whole set of sequential tests, or particular tests or series of tests, 
   can be triggered by issuing, in the ABINIT/tests directory, the command *./runtests.py*.
   (see the many capabilities of this scripts by issuing *./runtests.py --help*). 
   Other sets of calculations can be triggered, as described by issuing *make help*.
   The result will appear in a new subdirectory *Test_suite*

<!--
Additional note: 
The code uses a Pickle database (test_suite.cpkl) to store the set of objects representing the tests of the test suite. 
Remember to regenerate the database after any change to the TEST_INFO section or any modification of the configuration 
parameters of the test suite (__init__.py files). runtests.py provides the handy option -r (regenerate) 
to automatically regenerate the database before running the tests.
-->

Last but not least: are you sure that your modifications do not deteriorate the performance of the code 
in the regime where your modifications are not used?
You should inspect your modifications for both memory use and CPU time.

<!-- 
Include external files. Note that these files are also used by buildbot to 
provide hints when one of the tests fail so we have to keep them in separated files.
-->

{% include doc/developers/robodoc.doc.txt %}

{% include doc/developers/debug_make_parents %}

{% include doc/developers/debug_make_abiauty %}
