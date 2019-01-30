# WIP The new testsuite

## Presentation

### A name for the project ?

- ABINEAT: ABInit NEw physics Aware Testing
- NEAT: New Abinit Testing
- ABINUTS / NUTS: New Unit Testing System / New Unit Test Suite

### Principle

The goal of this project is to give to Abinit devellopers tools to create physics aware tests.
The mean to achieve this is a set of Fortran and Python tools to output structured data and
process it easily.

### Motivations

__Current situation__

In its current state the Abinit testing system is based on inputs files and associated reference output files.
A test consist in comparing the reference file with the outputfile of a fresh run in a pretty linear way,
only making difference between floating point numbers and other strings.

This approach is verry strict because it compare every single quantity in the main output file and report each time
a difference exceed a given tolerance. Moreover the line number (and therefor the number of iterations in every
loops printint in the output file) have to be exactly the same to allow the comparision to end properly and mark the
test as succeded.

Given the variaty of configurations of the test farm and the use of different parallelism configurations it may be
hard to create a stable test with decent tolerances and a lot of current test are really weak because of this limitation.

__Solution concepts__

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

### Proof of concept implementation

A few basic tools have already been implemented to check feasability.

On the Fortran side:
- a module based on C and Fortran code have been written for managing a map-like 
key-value pair list. Pair list have been chosen over hash table because in the practical case the O(n) cost of
finding a key is less problematic than the overhead in memory usage of hash tables because their will never 
be a lot of elements in the table.
The pair list can dynamicaly hold either real of integer values. It implement the basic getters and setters 
and an iterator system allowing to loop over the key-value pairs easily in Fortran.
- a module written in Fortran and using the pair list module have been written to easily output YAML into an
arbitrary file from Fortran code. This module is supposed to be the lowest level layer of the outputing system
and may be abstracted by higher layers. It provides basic routines for outputing YAML documents based on a
mapping at the root and containing well formated (as human readable as possible and valid YAML) 1D and 2D arrays
of numbers (integer or real), key-value mapping (using pair lists), list of key-value mapping and scalar fields.
Routines support tags for specifying special data structures.

On the python side:
- fldiff algorithm have been slightly modified to extract encountered YAML documents from the source and reserve
them for later treatment.
- Tools have been written to make easier the creation of new classes corresponding to YAML tags for adding new
features to the extracted datas. These tools are class decorators very easy to use.
- These tools have been used to create basic classes for futures tags, among other classes that directly convert
YAML list of numbers into NumPy arrays. These classes may be used as examples for the creation of furter tags.

