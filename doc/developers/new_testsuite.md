# WIP The new testsuite

## Presentation

### A name for the project ?

ABINEAT: ABInit NEw physics Aware Testing
NEAT: New Abinit Testing
NUTS: New Unit Testing System / New Unit Test Suite

### Principle

The goal of this project is to give to Abinit devellopers tools to create physics aware tests.
The mean to achieve this is a set of Fortran and Python tools to output structured data and
process it easily.

### Motivations

#### Current situation

In its current state the Abinit testing system is based on inputs files and associated reference output files.
A test consist in comparing the reference file with the outputfile of a fresh run in a pretty linear way,
only making difference between floating point numbers and other strings.

This approach is verry strict because it compare every single quantity in the main output file and report each time
a difference exceed a given tolerance. Moreover the line number (and therefor the number of iterations in every
loops printint in the output file) have to be exactly the same to allow the comparision to end properly and mark the
test as succeded.

Given the variaty of configurations of the test farm and the use of different parallelism configurations it may be
hard to create a stable test with decent tolerances and a lot of current test are really weak because of this limitation.

#### Solution concepts

The proposal of this project is to allow stricter tests on pertinent quantities, ignoring the rest. This can be achieved
by creating tests specialized for given physical quantities, taking in account there properties (Newton's Law,
energy conservation, well known equations...).

Of course all these tests would have to be designed nearly independantly for every interesting quantities and that is why
this project should focus on producing reusable tools to create these tests.

These tools would be used to create the tests associated with the more common quantities of ground state calculations (
composants of total energy, stress tensor and atomic forces). The test suite would be enriched later with the
participation of the Abinit community.

## Draft ideas

The outputing of data should use a machine readable format, preferably human readable but not necessarily.
This format could be YAML which is already used in a few situations in Abinit, JSON (for is great support)
or NetCDF (already used in Abipy).
Post processing tools may be available like plotting and computing facilities based on standard python
technologies like NumPy, SciPy, Matplotlib and Pandas.

For the Fortran side routines may be available to output given structures of data (matrices, scalars, arrays...) along with
meaningful metadata helping post processing and analysis.


