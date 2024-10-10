# Without even a single test, it is likely that an input variable would cease to work after a few years...

## So: how to add a new test in the test suite ?

The use of the test suite, from the user's point of view, has been introduced in the installation notes (see section 3 and 5). The following information complements the installation notes for the developer's point of view.

In order to introduce a test, one needs to:

- Provide **a new input file** in tests/v8/Input (or modify an existing one, possibly in another directory, but please do not suppress the existing capability testing!)
- Provide **a new reference file** in the corresponding tests/v8/Refs
- Insert the test in **the list of tests to be done**, by adding its name in tests/v8/__init__.py
- **Document the test** by adding a commented section inside the input file (edit an existing input file to follow its style)
- **Declare the pseudopotentials, the flow of actions, the files to be analyzed, the tolerances for the test**, inside this documentation.

For the tolerances, start with **tolnlines = 0, tolabs = 0.000e+00, tolrel = 0.000e+00**

The different fields control how strictly the test results will be analyzed:
- the maximum floating point difference found must be less than tolabs,
- the maximum relative difference less than tolrel, for each line individually,
- and tolnlines is the maximal number of lines found to be different between the output and reference file (within another tolerance that can be tuned by the option opt=-medium, opt=-easy or opt=-ridiculous, see the added section in other input files).

The scripts tests/Scripts/fldiff.pl and tests/Scripts/reportdiff.pl analyze the output.

If this procedure fails, contact the code maintainers in order adjust the test case. This might mean modifying the tolerances files. Unless you are really expert, let the maintainer do the final adjustment.

**Please, try to keep preferably the total time per test case to less than 10 seconds. 30 seconds is a maximum for most test cases. Going beyond this needs exceptional reasons.**

Supposing now that you have introduced a new test case. It can be used, with the other tests of the ABINIT test suite, through different channels:

- these tests can be triggered on the test farm (see the [buildbot 
  particulars](misc.md) Web page)
- locally (on your machine), the whole set of sequential tests, or particular tests or series of tests, can be triggered by issuing, in the ABINIT/tests directory, the command "./runtests.py" (see the many capabilities of this scripts by issuing ./runtests.py --help). Other sets of calculations can be triggered, as described by issuing "make help". The result will appear in a new subdirectory Test_suite

### Additional note

The code uses a Pickle database (test_suite.cpkl) to store the set of objects representing the tests of the test suite. Remember to regenerate the database after any change to the TEST_INFO section or any modification of the configuration parameters of the test suite (__init__.py files). runtests.py provides the handy option -r (regenerate) to automatically regenerate the database before running the tests.
