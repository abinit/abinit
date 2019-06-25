---
description: Design of the new Yaml-based test suite.
authors: TC, MG
---

# Yaml-based test suite: design principles and implementation details

This document discusses the new Yaml-based format used to represent physical results
in the main output file.
The connection with the ABINIT test suite and the syntax used to define tolerances, parameters and constraints in Yaml syntax
is also discussed.

## Presentation

### A name for the project

NEAT: NEw Abinit Test system

### Principle

The goal of this project is to provide ABINIT developers with reusable tools to implement physics-aware tests. 
With "physics-aware" we mean tests in which the developer can customize the tolerances and the logic used to 
compare numerical results thanks to the fact that the new parser
is aware of the context and of the meaning of the numerical values extracted from the output file.
The new infrastructure consists of a set of Fortran modules to output 
structured data in YAML format and Python code to parse the output files and analyze data.

### Motivations

__Situation in ABINITv8__

In ABINITv8 and previous versions, the ABINIT test suite is based on input files with the
corresponding reference output files.  Roughly speaking, an automatic test
consists in comparing the reference file with the output file in a line-oriented
fashion by computing differences between floating-point numbers without any
knowledge about the meaning and the importance of the numerical values.  This
approach is very rigid and rather limited because it obliges developers to use a
single (usually large) tolerance to account for possibly large fluctuations in
the intermediate results whereas the tolerance criteria should be ideally
applied only to the final results that are (hopefully) independent of details
such as hardware, optimization level and parallelism.

This limitation is clearly seen when comparing the results of iterative
algorithms.  The number of iterations required to converge, indeed, may depend
on several factors especially when the input parameters are far from convergence
or when different parallelization schemes or stochastic methods are used.  From
the point of view of code validation, what really matters is the final converged
value plus the time to solution if performance ends up being of concern.
Unfortunately, any line-by-line comparison algorithm will miserably fail in such
conditions because it will continue to insist on having the same number of
iterations (lines) in the two calculations to mark the test as succeeded.

It is evident that the approach used so far to validate new developments in
Abinit is not able to cope with the challenges posed by high-performance
computing and that smarter and more flexible approaches are needed to address
these limitations.

__Solution__

One of the goals of this project is to implement a python-based infrastructure
that allows developers to implement more rigorous tests on the portability of
important physical quantities while ignoring intermediate results, or using more relaxed
criteria for them. This goal is
achieved by providing a *declarative* interface that allows developers to define
the logic to be used to compare selected physical quantities.  We also provide a
declarative API to check that the results computed by the ab-initio code
satisfies fundamental physical and mathematical rules such as energy
conservation, Newton's third law, symmetry properties, etc.

That's not quite as simple as it sounds, especially because one of our goals is
to minimize the amount of (coding) work required to express such logic.  There
are therefore basic rules and design principles Abinit developers should be
aware of in order to take fully advantage of the new infrastructure.  In what
follows, we briefly describe the philosophy employed to implement the YAML-based
test procedure, the syntax used to compare outputs with references and the modifications required to
extend the framework.

<!--
Of course all these tests would have to be designed nearly independently for
every interesting quantities and that is why this project should focus on
producing reusable tools to create these tests.
These tools would be used to create the tests associated with the more common
quantities of ground state calculations (components of total energy, stress
tensor and atomic forces). The test suite would be enriched later with the
participation of the ABINIT community.
-->

!!! important

    The new Yaml-based testsuite relies on libraries that are not provided 
    by the python standard library. 
    To install these dependencies in *user* mode, use:


        pip install numpy pyyaml pandas --user

    If these dependencies are not available, the new test system will be
    disabled. Also you will have a warning saying so.


## Implementation details

The most important physical results are inserted in the main output file (aka
*ab_out*) inside machine-readable [YAML](https://en.wikipedia.org/wiki/YAML)
documents. A YAML document starts with three hyphens (---) followed by an
(optional) tag beginning with an exclamation mark (e.g. `!ETOT`).  Three periods
(...) signals the end of the document.  According to this syntax, one can easily
write a dictionary storing the different contributions to the total free energy
with:


```yaml
--- !ETOT
label               : Etot
comment             : Components of total free energy (in Hartree)
Kinetic energy      :  5.279019930263079807E+00
Hartree energy      :  8.846303409910728499E-01
XC energy           : -4.035286123400158687E+00
Ewald energy        : -1.095155892074578219E+01
PsPCore             :  1.549159282251551528E-01
Loc. psp. energy    : -3.356800625695788298E+00
NL psp energy       :  1.895544586353373973E+00
Etotal              : -1.012953488400904689E+01
Band energy         :  1.406569063148472631E+00
Total energy(eV)    : -2.756386620520307815E+02
...
```

Further details about the meaning of tags, labels and their connection with the
testing infrastructure will be given in the sections below.  For the time being,
it is sufficient to say that we opted for YAML because it is a *human-readable*
data serialization language already used in the log file to record important
events such as WARNINGs, ERRORs and COMMENTs (this is indeed the *protocol*
employed by AbiPy to monitor the status of Abinit calculations). 

Many programming languages, including python, provide support for reading and
writing YAML hence it is relatively easy to implement post processing tools
based on well-established python libraries for scientific computing such as
NumPy, SciPy, Matplotlib and Pandas.  Last but not least, writing YAML in
Fortran does not represent an insurmountable problem provided one keeps the
level of complexity of the YAML document at a *reasonable level*.  Our Fortran
implementation, indeed, supports only a subset of the YAML specifications:

* scalars
* arrays of one or two dimensions
* mappings containing scalars
* tables in CSV format (which is more like an extension of the standard)

!!! important

    YAML is not designed to handle large amount of data therefore it should
    not be used to represent large arrays for which performance is critical and human-readability 
    is lost by definition (do you really consider a YAML list with one thousand numbers human-readable?).
    Following this philosophy, YAML is supposed to be used to print the most important results 
    in the main output file and should not be considered as a replacement for binary netcdf files 
    when it comes to storing large data structures with lots of metadata.

    Note also that we do not plan to rewrite entirely the main output file in YAML syntax
    but we prefer to focus on those physical properties that will be used by the new test procedure
    to validate new developments.
    This approach, indeed, will facilitate the migration to the new YAML-based approach as only selected portions
    of the output file will be ported to the new format thus maintaining the look and the feel relatively 
    close to the previous *unstructured* format.

Several low-level tools have already been implemented during the prototyping process:

On the Fortran side, we provide:

- a Fortran module for managing dictionaries mapping strings to values. 
  Internally, the dictionary is implemented in C in terms of a list of key-value pairs.
  Lists have been chosen over hash tables because, in this particular context, 
  performance is not critical and the O(n) cost required to locate a key is negligible 
  given that there will never be a lot of elements in the table. 
  The pair list can dynamically hold either real, integer or string values. 
  It implements the basic getters and setters and an iterator system to allow looping 
  over the key-value pairs easily in Fortran.

- a Fortran module providing tools to manipulate variable length
  strings with a file-like interface (stream object)

- a Fortran module based on the two previous modules providing the API
  to easily output YAML documents from Fortran. 
  This module is supposed to be the lowest level layer of the output system and
  may be abstracted by higher layers. It provides basic routines to produce
  YAML documents based on a mapping at the root and containing well formatted
  (as human readable as possible and valid YAML) 1D and 2D arrays of numbers
  (integer or real), key-value mapping (using pair lists), list of key-value
  mapping and scalar fields.  Routines support tags for specifying special data
  structures.

- a higher level module called *m_neat* provides Fortran procedures to create specific
  documents associated to important physical properties. Currently routines
  for total energy components and ground state general results are implemented.

A more detailed example is provided [below](#creating-a-new-document-in-fortran).

<!--
For the Fortran side routines may be available to output given structures of
data (matrices, scalars, arrays...) along with meaningful metadata helping post
processing and analysis.
-->

As concerns the python implementation:

- The fldiff algorithm has been slightly modified to extract YAML
  documents from the output file and reserve them for later treatment.

- Tools have been written to facilitate the creation of new python classes
  corresponding to YAML tags for adding new features to the extracted data.
  These tools are very easy to use class decorators.

- These tools have been used to create basic classes for futures tags, among
  other classes that directly convert YAML list of numbers into NumPy arrays.
  These classes may be used as examples for the creation of further tags.

- A parser for test configuration have been added and all facilities to do
  tests are in place.

- A command line tool `testtools.py` to allow doing different manual actions
  (see Test CLI)

The new infrastructure has been designed with extensibility and ease-of-use in
mind. From the perspective of an Abinit developer, adding support for the
YAML-based approach requires two steps:

1. Implement the output of the YAML document in Fortran using the pre-existent
   API. Associate a label and possibly a tag to the new document.
2. Create a YAML configuration for the new document in tests that will produce
   it

If new tags are used (to handle new advanced structures) a third step is
required consisting in registering the tag and the associated class in the
python side. More about that can be found below.

An example will help clarify.  Let's assume we want to implement the output of
the `!ETOT` dictionary with the different components of the total free energy.
The workflow is split in two parts. The idea is to separate the computation
from the composition of the document, and separate the composition of the
document from the actual rendering.

In Fortran, we will have to write in the file of the computation:

```fortran
!! Header of the module
use m_neat, only: neat_etot
use m_pair_list, only pair_list
...

!! header of the routine
real(kind=dp) :: etot
type(pair_list) :: e_components
...

!! body of the routine
call e_components%set("comment", s="Total energie and its components")
call e_components%set("Total Energy", r=etot)
...

!! footer of the routine, once all data have been stored
call neat_etot(e_components)
```

In *47_neat/m\_neat.F90* we will implement *neat\_etot*:
```fortran
subroutine neat_etot(components, abiout)
type(pair_list), intent(in) :: components
integer, intent(in) :: abiout
    call yaml_single_dict("Etot", "", components, 35, 500, tag="ETOT", width=20, file=abiout, real_fmt='(ES20.13)')
    !! 35 -> max size of the label, needed to extract from the pairlist, 500 -> max size of the strings, needed to extract the comment
    !! width -> width of the field name side, permit a nice alignment of the values
end subroutine
```

## Test specification draft

### Declaration from ABINIT input file

In the standard approach, the parameters governing the execution of the test
are specified in the `TEST_INFO` section located at the end of the file.
The options are given in the [INI file format](https://en.wikipedia.org/wiki/INI_file).
The integration of the new YAML-based tests with the pre-existent infrastructure 
is obtained via two modifications of the current specifications.
More specifically:

(1) the *files_to_test* section now accepts the optional argument *use_yaml*. 
   Possible values are:

   - "yes" --> activate YAML based test
   - "no" -->  do not use YAML test (default)
   - "only" --> use YAML based test but not legacy fldiff algorithm

(2) a new (optional) section `[yaml_test]` is added with two possible fields: 

   - *test* --> the value is used as the YAML test specification source, it may
   be used for really short configuration heavily relying on the default.
   - *file* --> the value is the path to the file containing the YAML test
   specification. The path is relative to the Abinit input file. A natural choice
   would be to use the path "./t21.yaml" associated to the input file "t21.in"
   unless we decide to create a dedicated directory.

!!! important
    Until the basic document list is considered stable enough the printing of
    YAML documents is disabled by default. To enable it, add `use_yaml 1` to your
    input file (in the normal Abinit input variables part).

### Test specification syntax

To explain how to build a test specification I will go step-by-step through
the definition of the _paral[86]_ test.

First we have to introduce the main concepts.

Data tree
: The processed data is a YAML document. Therefore it can be seen as a tree
  whose nodes have a label and leaf are scalars or special data structures
  identified by a tag (not all tags mark a leaf). The top-level nodes are the
  YAML documents themselves.

Config tree
: The configuration also takes the form of a tree. Its nodes are
  _specializations_ and its leaf are _parameters_ or _constraints_.
  Its structure matches the structure of the _data tree_ so one can define rules
  (constraint and parameters) that apply to a specific place of the _data tree_.

Specialization
: The rules defined under a specialization will apply only on the matching node
  of the _data tree_ and its children.

Constraint
: A constraint is a condition one impose for the test to succeed. Constraints
  can apply to leafs of the data tree or to nodes depending of the nature of the
  constraint.

Parameter
: A parameter is a value that can be used by constraints to modify their
  behavior.

Iteration state
: An iteration state describe how much iterations of each possible level have
  happened in the run (ex: idtset=2, itimimage=not used, image=5, time=not used).
  It give information on the current state of the run. Documents are implicitly
  associated to their iteration state. These informations are available to
  the test engine through special YAML documents using the `IterStart` tag.

!!! tip

    To get the list of constraints and parameters, run:
    
        ~abinit/tests/testtools.py explore
    
    and type `show *`. You can then type for example `show tol_eq` to learn more
    about a specific constraint or parameter.

A few conventions on documents writing
: The label field appear in all data document. It should be a unique identifier
  at the scale of an _iteration state_. The tag is not necessarily unique, it 
  describe the structure of the document and there is no need to use it unless
  special logic have to be implemented for the document. The comment field is
  optional but is appreciated where the purpose of the document is not obvious

Let's start with a minimalistic example in which we compare the components of
the total free energy.

We will start by comparing the values of each components with both absolute
tolerance of 1.0e-7 Ha.

The document holding this piece of data have the _label_ `Etot` so the YAML file
will look like this:

```yaml
Etot:
    tol_abs: 1.0e-7
```

The `tol_abs` keyword is a _constraint_. It will check for every child of
`Etot` that the value does not differ by more than 1.0e-7 Ha with the reference. 
If one does, the test will fail and the report will tell you which component is wrong.

In the `Etot` document one can find the total energy in eV. Since the unit is
different the absolute tolerance have not the same impact on the precision.
We want to achieve the same relative precision on this term but we cannot
achieve the same absolute precision. Though we will __specialize__ the rules
and use the `tol_rel` constraint. The name of the field holding the total energy
in eV in the `Etot` document is `Total energy (eV)`. We could do something like
this:

```yaml
Etot:
    tol_abs: 1.0e-7
    Total energy (eV):
        tol_rel: 1.0e-10
```

However the `tol_abs` constraint defined in `Etot` is _inherited_ by `Total
energy (eV)` which mean we will not only check `tol_rel` with a relative tolerance of
1.0e-10 but also `tol_abs` with an absolute tolerance of 1.0e-7. Most of the time it is
what we need, even though here, it is not. So we will just use a different
tolerance for `tol_abs` in `Total energy (eV)`.

```yaml
Etot:
    tol_abs: 1.0e-7
    Total energy (eV):
        tol_abs: 1.0e-5
        tol_rel: 1.0e-10
```

Now we achieve the same relative precision and the test does not fail because
of the looser absolute precision of the total energy in eV.

The `Etot` document is the simplest possible document. It only contains fields
with real values. We now will have a look at the `results_gs` document. It comes
from the structure of the same name in ABINIT code. In the ABINIT output file it
look like this:

```yaml
---
label     : results_gs
comment   : Summary of ground states results.
natom     :        5
nsppol    :        1
cut       : {"ecut":   1.20000000000000000E+01, "pawecutdg":   2.00000000000000000E+01, }
convergence: {
    "deltae":   2.37409381043107715E-09, "res2":   1.41518780109792898E-08, 
    "residm":   2.60254842131463755E-07, "diffor": 0.00000000000000000E+00, 
}
etotal    :  -1.51507711707660150E+02
entropy   :   0.00000000000000000E+00
fermie    :   3.09658145725792422E-01
stress tensor: !Tensor
- [  3.56483996349480498E-03,   0.00000000000000000E+00,   0.00000000000000000E+00, ]
- [  0.00000000000000000E+00,   3.56483996349480151E-03,   0.00000000000000000E+00, ]
- [  0.00000000000000000E+00,   0.00000000000000000E+00,   3.56483996349478416E-03, ]

cartesian forces: !CartForces
- [ -0.00000000000000000E+00,  -0.00000000000000000E+00,  -0.00000000000000000E+00, ]
- [ -0.00000000000000000E+00,  -0.00000000000000000E+00,  -0.00000000000000000E+00, ]
- [ -0.00000000000000000E+00,  -0.00000000000000000E+00,  -0.00000000000000000E+00, ]
- [ -0.00000000000000000E+00,  -0.00000000000000000E+00,  -0.00000000000000000E+00, ]
- [ -0.00000000000000000E+00,  -0.00000000000000000E+00,  -0.00000000000000000E+00, ]
...
```

Here we have a various examples of what can be found in a document. We have real
and integer fields, mappings/dictionaries and 2D arrays.

To check this document in our test we have to add a specialization for it along
side the `Etot` one. It will look like this:

```yaml
Etot:
    tol_abs: 1.0e-7
    Total energy (eV):
        tol_abs: 1.0e-5
        tol_rel: 1.0e-10
results_gs:
    tol_rel: 1.0e-8
```

For simplicity sake I will only write the `results_gs` part in next examples.
First we want to check the integers and real values and so we use `tol_rel`.
However the content of the `convergence` field represent residues. It does not
make sense to check those against the reference. What we really want is make
sure that they are below a given ceil. This is what the `ceil` constraint is
for:

```yaml
results_gs:
    tol_rel: 1.0e-8
    convergence:
        ceil: 3.0e-7
```

Now the test will fail if one of the components of the `convergence` mapping is
above 3.0e-7. The `ceil` constraint have the particularity to disable `tol_rel`
and `tol_abs` defined before because they are mutually exclusive.

!!! tip

    Within the explore shell `show ceil` will tell you, among other informations,
    what are the constraints disabled by the use of `ceil` in the _exclude_
    field.

Fields with the `!Tensor` tags are leafs of the tree. Then the tester routine
won't try to compare each individual coefficient with `tol_rel`. However we
still want to check that it does not change too much. For that purpose we use the
`tol_vec` constraint which apply to all arrays derived from `BaseArray` (most
arrays with a tag). `BaseArray` let us use the capabilities of Numpy arrays with
YAML defined arrays. `tol_vec` check the euclidian distance between reference
and tested array. Since we also want to apply this constraint to
`cartesian_force` we will define the constraint at the top level of `results_gs`.

```yaml
results_gs:
    tol_rel: 1.0e-8
    convergence:
        ceil: 3.0e-7
    tol_vec: 1.0e-5
```

We now have a basic test that check the values that matter with a fair
precision.
However this isn't satisfying yet. Indeed the test _paral[86]_ like a lot
of Abinit tests use two datasets to perform two very different computations.
The first dataset is dedicated to the computation of a good density in DFT+LDA. 
The second dataset use this density to start a computation with DMFT. Since the
purpose of the two dataset are very different, different convergence criterion
are defined in the input file. Then we would like to create different configurations
 for each dataset. This is possible thanks to __filters__.

A filter is a way to associate a specific configuration to a set of _iteration
state_. A filter is defined in a separated section of the test configuration
under the node `filters`. We will add the following to our config file to
declare two filters:

```yaml
filters:
    lda:
        dtset: 1
    dmft:
        dtset: 2
```

Here we simply said that we associate the label `lda` to a filter matching all
documents created in the first dataset and we associate the label `dmft` to a
filter matching all document created in the second dataset. This is the simplest
filter declaration possible. More on filters declarations in
[here](yaml_tools_api#filters-api).

Now we have to use our filters. First of all we will associate the configuration
we already wrote to the `lda` filter so we can have a different configuration
for the second dataset. The YAML file will now became

```yaml
lda:
    Etot:
        tol_abs: 1.0e-7
        Total energy (eV):
            tol_abs: 1.0e-5
            tol_rel: 1.0e-10

    results_gs:
        tol_rel: 1.0e-8
        convergence:
            ceil: 3.0e-7
        tol_vec: 1.0e-5

filters:
    lda:
        dtset: 1
    dmft:
        dtset: 2
```

By putting the configuration under the `lda` node we state that those rules only
apply to the first dataset. We will then create a new `dmft` node and create a
configuration following the same procedure than before. We end up with something
like this:

```yaml
lda:
    Etot:
        tol_abs: 1.0e-7
        Total energy (eV):
            tol_abs: 1.0e-5
            tol_rel: 1.0e-10

    results_gs:
        tol_rel: 1.0e-8
        convergence:
            ceil: 3.0e-7
        tol_vec: 1.0e-5

dmft:
    tol_abs: 2.0e-8
    tol_rel: 5.0e-9

    results_gs:
        convergence:
            ceil: 1.0e-6
            diffor:
                ignore: true
        fermie:
            tol_abs: 1.0e-7
            tol_rel: 1.0e-8

        stress tensor:
            ignore: true
    Etot:
        Total energy eV:
            tol_abs: 1.0e-5
            tol_rel: 1.0e-8

    Etot DC:
        tol_abs: 1.0e-7
        Total DC energy eV:
            tol_abs: 1.0e-5
            tol_rel: 1.0e-8

filters:
    lda:
        dtset: 1
    dmft:
        dtset: 2
```

## Command line interface

The `~abinit/tests/testtools.py` script provides a command line interface to facilitate the creation of new tests. 
The script uses the syntax:

    ./testtools.py COMMAND [options]

Run the script without arguments to get the list of possible commands and use

    ./testtools.py COMMAND --help

to display the options supported by `COMMAND`.

The available commands are:

fldiff

:   Command line interface to the *fldiff.py* module. 
    It may be useful to compare output and reference files without running ABINIT each time 
    like *runtests.py* would do.
    It allows the user to provide the same parameters that are passed by the
    *testsuite.py* when *runtests.py* is used.


explore

:   This tool allows the user to *explore* and validate a test configuration
    file. It provides a shell like interface in which the user can move around the
    tree of the configuration file and print the constraints defined by the test. It
    also provides documentation about constraints and parameters via the *show*
    command. 

## Extending the test suite

### Entry points

There are three main entry points of growing complexity in this system.

The first one is the yaml configuration file intended to all Abinit developers. It
does not require any Python knowledge, only a basic comprehension of the conventions
used for writing tests. Being fully *declarative* (no logic) it should be quite easy
to learn its usage from the available examples.

The second one is the *~abinit/tests/pymods/yaml_tests/conf_parser.py* file. It
contains declarations of the available constraints and parameters. A basic
python understanding is required in order to modify this file. Comments and doc strings
should help users to grasp the meaning of this file. More details on 
[here](./yaml_tools_api#constraints-and-parameters-registration)

The third is the file *~abinit/tests/pymods/yaml_tests/structures/*. It
defines the structures used by the YAML parser when encountering a tag (starting
with !), or in some cases when reaching a given pattern (__undef__ for example).
The *structures* directory is a package organised by features (ex: there is a
file *structures/ground_state.py*). Each file define structures for a given
feature. All files are imported by the main script *structures/__init__.py*.
Even if the abstraction layer on top of the _yaml_ module should help, it is
better to have a good understanding of more "advanced" python concepts like
_inheritance_, _decorators_, _classmethod_ etc.

### The special constraints equation, equations and callback

`equation`, `equations` and `callback` are special constraints because their
actual effects are defined directly in the configuration file. They are made to
provide a bit more flexibility to the configuration file without diving into
python code.

equation and equations

: These two are sisters. `equation` take a string as a value. This string will
  be interpreted as a python expression that must result in a number. The
  absolute value of this number will be compared to the value of the `tol_eq`
  parameter and if `tol_eq` is greater the test will succeed. The expression can
  also result in a Numpy array. In this case it is the euclidean norm of the
  array that will be compared to `tol_eq` value. `equations` works exactly the
  same but has a list of string as value. Each string is a different expression
  that will be tested independently from the others. In both case the tested
  object can be referred as `this` and the reference object can be referred as
  `ref`.

Dumb example:
```yaml
Etot:
    tol_eq: 1.0e-6
    equation: 'this["Etotal"] - this["Total energy(eV)"]/27.2114'
```

callback

: This one require a bit of python coding since it will use a method of the
  structure it is defined in. Suppose we have a tag `!AtomSpeeds` associated
  to a class `AtomSpeeds` and to a document labeled `Atomic speeds` in the data
  tree. The `AtomSpeeds` class have a method `not_going_anywhere` that checks
  that the atoms are not going to try to leave the box. We would like to
  pass some kind of tolerance `d_min` the minimal distance atoms can approach
  the border of the box. The signature of the method have to be
  `not_going_anywhere(self, tested, d_min=DEFAULT_VALUE)` and should return
  `True`, `False` or an instance of `FailDetail` (see __Add a new constraint__
  for explanations about those). Note that `self` will be the reference
  instance. We can then use it by with the following configuration:

```yaml
Atomic speeds:
    callback:
        method: not_going_anywhere
        d_min: 1.0e-2
```

### Add a new tag

Pyyaml offer the possibility to directly convert some YAML structures to a
Python class using a tag. To register a new tag one should edit the file
`~abinit/tests/pymods/yaml_tools/structures.py`. In this file are defined several
classes that are decorated with one of `@yaml_map`, `@yaml_scalar`, `@yaml_seq`
or `@yaml_auto_map`. Those decorator are the functions that actually register
the class as a known tag. Whether you should use one or another depend on the
structure of the data in YAML. Is it a mapping/dictionary, a scalar or a
list/sequence ?

In most cases one wants `tester` to browse the children of the strucure and apply
relevant constraints on it and to access attributes through their original name
from the data. In this case the simpler way to register a new tag is to use
`yaml_auto_map`. For example to register a tag ETOT that simply register all
fields from the data tree and let the tester check them with `tol_abs` or
`tol_rel` we would put the following in
*~abinit/tests/pymods/yaml_tools/structures/ground_state.py*:

```python
@yaml_auto_map
class Etot(object):
    __yaml_tag = 'ETOT'
```

ETOT does not match the python naming conventions so we use `Etot` as a class
name, however since we want to use the tag `!ETOT` we use the special attribute
`__yaml_tag` to declare our custom tag. If this attribute is not used the name
of the class is the default tag.

`yaml_auto_map` does several things for us:
- it gives the class a `dict` like interface by defining relevant methods, which
  allows tester to browse children
- it registers the tag in YAML parser
- it automatically register all attributes found in the data tree as attributes
  of the class instance. These attributes are accessible through the attribute
  syntax (ex: `my_object.my_attribute_unit`) with a normalized name (basically remove
  characters that cannot be in a python identifier like spaces and ponctuation)
  or through the dictionary syntax with their original name (ex: `my_object['my
attribute (unit)']`)

Sometimes one wants more control over the building of the class instance. This
is what `yaml_map` is for.
Let suppose we still want to register tag but we want to select only a subset of
the components for example. We will use `yaml_map` to gives us control over the
building of the instance. This is done by implementing the `from_map` class
method (a class method is a method that is called from the class instead of
being called from an instance). This method take in argument a dictionary built
from the data tree by the YAML parser and should return an instance of our
class.

```python
@yaml_map
class Etot(object):
    __yaml_tag = 'ETOT'
    def __init__(self, kin, hart, xc, ew):
        self.kinetic = kin
        self.hartree = hart
        self.xc = xc
        self.ewald = ew

    @classmethod
    def from_map(cls, d):
        # cls is the class (Etot here but it can be
        # something else if we subclass Etot)
        # d is a dictionary built from the data tree
        kin = d['Kinetic energy']
        hart = d['Hartree energy']
        xc = d['XC energy']
        ew = d['Ewald energy']
        return cls(kin, hart, xc, ew)
```

Now we fully control the building of the structure, however we lost the ability
for the tester to browse the components to check them. If we want to only make
our custom check it is fine. For example we can define a method
`check_components` and use the `callback` constraint like this:

In `~abinit/tests/pymods/yaml_tools/structure/ground_state.py`
```python
@yaml_map
class Etot(object):
    __yaml_tag = 'ETOT'

    ...  # same code as the previous example

    def check_components(self, tested):
        # self will always be the reference object
        return (
            abs(self.kinetic - other.kinetic) < 1.0e-10
            and abs(self.hartree - other.hartree) < 1.0e-10
            and abs(self.xc - other.xc) < 1.0e-10
            and abs(self.ewald - other.ewald) < 1.0e-10
        )
```

In the configuration file
```
Etotal:
    callback:
        method: check_components
```

However it can be better to give back the control to `tester`. For that purpose
we can implement the method `get_children` that should return a dictionary of the
data to be checked automatically. `tester` will detect it and use it to check
what you gave him.

```python
@yaml_map
class Etot(object):
    __yaml_tag = 'ETOT'

    ...  # same code as the previous example

    def get_children(self):
        return {
            'kin': self.kinetic,
            'hart': self.hartree,
            'xc': self.xc,
            'ew': self.ewald
        }
```

Now `tester` will be able to apply `tol_abs` and friends to the components we
game him.

If the class has a __complete__ dict-like read interface (`__iter__` yielding
keys, `__contains__`, `__getitem__`, `keys` and `items`) then it can have a
class attribute `is_dict_like` set to `True` and it will be treated as any other
node (it not longer need `get_children`). `yaml_auto_map` registered classes
automatically address these requirements.

`yaml_seq` is analogous to `yaml_map` however `to_map` became `to_seq` and the
YAML source data have to match the YAML sequence structure (either `[a, b, c]`
or
```yaml
- a
- b
- c
```
)

The argument passed to `to_seq` is a list. If one wants `tester` to browse the
elements of the resulting object one can either implement a `get_children`
method or implement the `__iter__` python special method.

If for some reason a class have the `__iter__` method implemented but one does
__not__ want tester to browse its children (`BaseArray` and its subclasses are
such a case) one can defined the `has_no_child` attribute and set it to `True`.
Then tester won't try to browse it. Strings are a particular case. They have the
`__iter__` method but will never been browsed.

To associate a tag to anything else than a YAML mapping or sequence one can use
`yaml_scalar`. This decorator expect the class to have a method `from_scalar`
that takes the raw source as a string in argument and return an instance of the
class. It can be used to create custom parsers of new number representation.
For example to create a 3D vector with unit tag:
```python
@yaml_scalar
class Vec3Unit(object):
    def __init__(self, x, y, z, unit):
        self.x, self.y, self.z = x, y, z
        self.unit = unit
    
    @classmethod
    def from_scalar(cls, raw):
        sx, sy, sz, unit, *_ = raw.split()  # split on blanks
        return cls(float(sx), float(sy), float(sz), unit)
```

With that new tag registered the YAML parser will happily parse something like
```yaml
kpt: !Vec3Unit 0.5 0.5 0.5 Bohr^-1
```

Finally when the scalar have a easily detectable form one can create an implicit
scalar. An implicit scalar have the advantage to be detected by the YAML parser
without the need of writing the tag before as soon as it match a given regular
expression. For example to parse directly complex numbers we could (naively) do
the following:

```yaml
@yaml_implicit_scalar
class YAMLComplex(complex):
    # regular expression matching a (quite rigid) writing of complex number
    yaml_pattern = r'[+-]?\d+\.\d+) [+-] (\d+\.\d+)i'

    # this have nothing to do with yam_implicit_scalar, it is just the way to
    # create a subclass of a python *native* object
    @staticmethod
    def __new__(*args, **kwargs):
        return complex.__new__(*args, **kwargs)

    @classmethod
    def from_scalar(cls, scal):
        # few adjustment to match the python restrictions on complex number
        # writing
        return cls(scal.replace('i', 'j').replace(' ', ''))
```

### Creating a new document in Fortran

Developers are expected to browse the sources of `m_neat` and `m_yaml_out` to
have a comprehensive overview of available tools. Indeed the routines are
documented and commented and creating a hand-written reference is likely to go
out of sync quicker than on-site documentation. Here we show a little example to
give the feeling of the process.

The simplest way to create a YAML document have been introduced above with the
use of `yaml_single_dict`. However this method is limited to scalars only. It is
possible to put 1D and 2D arrays, dictionaries and even column based data.

The best way to create a new document is to have a routine
`neat_my_new_document` that will take all required data in argument and will do
all the formating work at once.

This routine will declare a `stream_string` object to build the YAML document
inside, open the document with `yaml_open_doc`, fill it, close it with
`yaml_close_doc` and finally print it to the output file with `wrtout_stream`.

Here come a basic skeleton of `neat` routine:

```fortran
subroutine neat_my_new_document(data_1, data_2,... , iout)
  !! declare your pieces of data... (arrays, numbers, pair_list...)
  ...
  integer,intent(in) :: iout  !! this is the output file descriptor
!Local variables-------------------------------
  type(stream_string) :: stream

  !! open the document
  call yaml_open_doc('my label', 'some comments on the document', stream=stream)
  
  !! fill the document
  ...
  
  !! close and output the document
  call yaml_close_doc(stream=stream)
  call wrtout_stream(stream, iout)
end subroutine neat_my_new_document
```

Suppose we want a 2D matrix of real number with the tag `!NiceMatrix` in our document for the field name 'that matrix'we will add the following to the middle section:
```fortran
call yaml_add_real2d('that matrix', dimension_1, dimension_2, mat_data, tag='NiceMatrix, stream=stream)
```

Other `m_yaml_out` routines provide similar interface: first the label, then the
data and its structural metadata, then a bunch of optional arguments for
formating, adding tags, tweaking spacing etc...


## Coding rules

This section discusses the basic rules that should be followed when writing Yaml documents in Fortran.
In a nutshell:

* no savage usage of write statements!
* In the majority of the cases, Yaml documents should be opened and closed in the same routine!
* A document with a tag is considered a standardized document i.e. a document for which there's 
  an official commitment from the ABINIT community to maintain backward compatibility.
  Official Yaml documents may be used by third-party software to implement post-processing tools.
* To be discussed: `crystal%yaml_write(stream, indent=4)`


## Further developments

This new Yaml-based infrastructure can be used as building block to implement the 
high-level logic required by more advanced integration tests such as:

Parametrized tests

: Tests in which multiple parameters are changed either at the level of the input variables 
  or at the MPI/OpenMP level.
  Typical example: running calculations with [[useylm]] in [0, 1] or [[paral_kgb]] = 1 runs 
  with multiple configurations of [[npfft]], [[npband]], [[npkpt]].

Benchmarks

: Tests to monitor the scalability of the code and make sure that serious bottlenecks are not 
  introduced in trunk/develop when important dimension are increased (e.g.
  [[chksymbreak]]  > 0 with [[ngkpt]] > 30\*\*3).
  Other possible applications: monitor the memory allocated to detect possible regressions.

Interface with AbiPy and Abiflows

: AbiPy has its own set of integration tests but here we mainly focus on the python layer
  without testing for numerical values. 
  Still it would be nice to check for numerical reproducibility, especially when it comes to 
  workflows that are already used in production for high-throughput applications (e.g. DFPT). 
