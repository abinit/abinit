# WIP The new testsuite

## Presentation

### A name for the project ?

- ABINEAT: ABInit NEw physics Aware Testing
- NEAT: New Abinit Testing
- ABINUTS / NUTS: New Unit Testing System / New Unit Test Suite

### Principle

The goal of this project is to give to Abinit devellopers tools to create
physics aware tests.  The mean to achieve this is a set of Fortran and Python
tools to output structured data and process it easily.

### Motivations

__Current situation__

In its current state the Abinit testing system is based on inputs files and
associated reference output files.  A test consist in comparing the reference
file with the outputfile of a fresh run in a pretty linear way, only making
difference between floating point numbers and other strings.

This approach is verry strict because it compare every single quantity in the
main output file and report each time a difference exceed a given tolerance.
Moreover the line number (and therefor the number of iterations in every loops
printint in the output file) have to be exactly the same to allow the
comparision to end properly and mark the test as succeded.

Given the variaty of configurations of the test farm and the use of different
parallelism configurations it may be hard to create a stable test with decent
tolerances and a lot of current test are really weak because of this limitation.

__Solution concepts__

The proposal of this project is to allow stricter tests on pertinent quantities,
ignoring the rest. This can be achieved by creating tests specialized for given
physical quantities, taking in account there properties (Newton's Law, energy
conservation, well known equations...).

Of course all these tests would have to be designed nearly independantly for
every interesting quantities and that is why this project should focus on
producing reusable tools to create these tests.

These tools would be used to create the tests associated with the more common
quantities of ground state calculations ( composants of total energy, stress
tensor and atomic forces). The test suite would be enriched later with the
participation of the Abinit community.


## Test specification draft

### Declaration from abinit input file

Abinit input files for tests already embed informations for the test suite in a
specific TEST\_INFO section at the end of the file. This section is structured
with a INI like syntax (the syntax define by the configparser module of
Python).  The introduction of YAML based test will be done with two
modifications to the current test specification:

- the *files\_to\_test* option will hold a new optional argument *use_yaml*
  that will hold one of "yes" (default, use YAML based test), "no" (do not use
  YAML test) or "only" (use YAML based test but not legacy fldiff algorithm)
- a new section *[yaml_test]* will optionally be added with two options: *test*
  whom value is an embed test specification and *file* which would be a path to
  the file containing the test specification

### Test specification syntax

The test specification will also use YAML syntax. It use three basic concepts:

- _constraint_ is an actual check of the data
- _parameter_ is an arbitrary value passed to some constraint and that can
  modify there behaviour
- _specialization_ is a way to override already defines _constraints_ and
  _parameters_ for a specific node of the data tree (and its children)

The test config file consists of a single YAML document structured as a
mapping. Each element of the mapping can be a constraint, a parameter of a
specialization. The distinction is made by looking at the list of known
parameters and constraints. The value of a parameter or a constraint depend on
its definition. The value of a specialization is itself a mapping filled with
parameters, constraints and specialization.

The first level of specialisation match the documents label value from tested files.
The next levels match the different attributes of the document.

To get the list of constraints and parameters run the program
`~abinit/tests/testcli.py shell` and type `show *`. You can then type for
example `show tol_eq` to learn more about a specific constraints or parameter.

Constraints and parameters have several properties that define their behaviour.
Most of the time constraints and parameters and apply to all children level,
they are inherited, though some of them apply only at the level where they are
define. Constraints can have a list of parameters they use. If these parameters
are not available at the current level (they are not defined in any higher
level or they are not herited) the default value of the parameter is used.
Constraints can be mutually exclusive which mean that when a constraint is defined,
all constraints define in higher level this one excludes are hidden to the
current level and its children.

Here is an exemple of a possible file:
```yaml
Etot:  # this define rules for the document Etot
    tol: 1.0e-3  # this say that for all numerical values encountered in
                 # the document Etot we will check a tolerance of 1.0e-3
                 # both on absolute an relative differences
    Kinetic energy:  # this override the tol constraints for the Kinetic enrgy
                     # attribute
        tol: 1.0e-10

results_gs:
    tol: 1.0e-5
    convergence:
        ceil: 5.0e-8  # here we define a ceil which exclude tol
                      # instead of a test between ref and tested file
                      # we will just check that values in convergence
                      # are below 5.0e-8
        residm:
            ceil: 5.0e-6

    entropy:
        ceil: 1.0e-10

    stress tensor:
        tol: 1.0e-5
        tensor_is_symetric: True  # check that tensor is symetric

    cartesian_force:
        tol_eq: 1.0e-8  # tol_eq is a parameter used by equation
        equation: "this.sum(axis=0)"  # 2-norm of sum of forces is under 1e-8

    etotal:
        tol: 1.0e-8
```

## Proof of concept implementation

The outputing of data should use a machine readable format, preferably human
readable but not necessarily.  This format could be YAML which is already used
in a few situations in Abinit, JSON (for its great support) or NetCDF (already
used in Abipy).  Post processing tools may be available like plotting and
computing facilities based on standard python technologies like NumPy, SciPy,
Matplotlib and Pandas. Currently the prototype integrate YAML symtax and NumPy
arrays.

For the Fortran side routines may be available to output given structures of
data (matrices, scalars, arrays...) along with meaningful metadata helping post
processing and analysis.

A few basic tools have already been implemented to check feasability.

On the Fortran side:

- a module based on C and Fortran code have been written for managing a map-like
  key-value pair list. Pair list have been chosen over hash table because in the
  practical case the O(n) cost of finding a key is less problematic than the
  overhead in memory usage of hash tables because their will never be a lot of
  elements in the table.  The pair list can dynamicaly hold either real of
  integer values. It implement the basic getters and setters and an iterator
  system allowing to loop over the key-value pairs easily in Fortran.
- a module written in Fortran provides tools to manipulate variable length strings 
  with a file like interface.
- a module written in Fortran and using the two previous module have been written
  to easily output YAML into an arbitrary file from Fortran code. This module is
  supposed to be the lowest level layer of the outputing system and may be
  abstracted by higher layers. It provides basic routines for outputing YAML
  documents based on a mapping at the root and containing well formated (as
  human readable as possible and valid YAML) 1D and 2D arrays of numbers
  (integer or real), key-value mapping (using pair lists), list of key-value
  mapping and scalar fields.  Routines support tags for specifying special data
  structures.
- a higher level module called m\_neat provide routines for creating specific
  documents associated with interesting physical properties. Currently routines
  for total energy components and ground state general results are implemented.

On the python side:

- fldiff algorithm have been slightly modified to extract encountered YAML
  documents from the source and reserve them for later treatment.
- Tools have been written to make easier the creation of new classes
  corresponding to YAML tags for adding new features to the extracted datas.
  These tools are very easy to use class decorators.
- These tools have been used to create basic classes for futures tags, among
  other classes that directly convert YAML list of numbers into NumPy arrays.
  These classes may be used as examples for the creation of furter tags.
- A parser for test configuration have been added and all facilities to do tests
  are in place.
- A command line tool testcli.py allow to do different manual actions (see Test CLI)

## Test CLI

A command line tool `~abinit/tests/testcli.py` provide tools to work on writing tests.
It provide help if run without argument.
The available sub commands are described here.

### Diff

The __diff__ subcommand provide a command line interface to the fldiff.py
module. It may be usefull to compare output and reference file without running
Abinit each time like runtest.py would do.

It allow the user to provide the same parameters that are passed by the
testsuite.py when runtest.py is used.

Use `~abinit/tests/testcli.py diff --help` for more informations.

### Shell

This tool provide is helpful to explore a test configuration file. It provide a
shell like interface where the user can move around the tree of the
configuration file, seeing what constraints define the test.  It also provide
documentation about contstraints.
Run `~abinit/tests/testcli.py shell` to use it.
