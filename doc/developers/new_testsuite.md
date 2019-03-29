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

### A name for the project ?

- ABINEAT: ABInit NEw physics Aware Testing
- NEAT: New Abinit Testing
- ABINUTS / NUTS: New Unit Testing System / New Unit Test Suite

### Principle

The goal of this project is to provide ABINIT developers with reusable tools to implement physics-aware tests. 
With "physics-aware" we mean tests in which the developer can customize the tolerances and the logic used to 
compare numerical results thanks to the fact that the new parser
is aware of the context and of the meaning of the numerical values extracted from the output file.
The new infrastructure consists of a set of Fortran modules to output 
structured data in YAML format and Python code to parse the output files and analyze data.

### Motivations

__Current situation__

In its current state, the ABINIT test suite is based on input files with the
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
important physical quantities while ignoring intermediate results.  This goal is
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
test suite, the syntax used to define tests and the modifications required to
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


## Implementation details

The most important physical results are written in the main output file (aka
*ab_out*) inside machine-readable [YAML](https://en.wikipedia.org/wiki/YAML)
documents.  A YAML document starts with three hyphens (---) followed by an
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
* mapping containing scalars
* tables in CSV format

A more detailed discussion about the Fortran API is given in the XXX section.

!!! important

    YAML is not designed to handle large amount of data therefore it should
    not be used to represent large arrays for which performance is critical and human-readability 
    is lost by definition (do you really consider a YAML list with one thousand numbers human-readable?).
    Following this philosophy, YAML is supposed to be used to print the most important results 
    in the main output file and should not be considered as a replacement for binary netcdf files 
    when it comes to storing large data structures with lots of metadata.

    Note also that we do not plan to rewrite entirely the main output file in YAML syntax
    but we prefer to focus on those physical properties that will be used by the new test suite 
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
  strings with a file like interface (stream object)

- a Fortran module based on the two previous modules providing the API
  to easily output YAML documents from Fortran. 
  This module is supposed to be the lowest level layer of the output system and
  may be abstracted by higher layers. It provides basic routines to produce
  YAML documents based on a mapping at the root and containing well formatted
  (as human readable as possible and valid YAML) 1D and 2D arrays of numbers
  (integer or real), key-value mapping (using pair lists), list of key-value
  mapping and scalar fields.  Routines support tags for specifying special data
  structures.

- a higher level module called *m_neat* provides Fortran procedures to creat specific
  documents associated to important physical properties. Currently routines
  for total energy components and ground state general results are implemented.

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

==TODO== Add Pseudo code with examples showing how to write ETOT in Fortran, how
to define and register a python YAML object associated to the ETOT tag.  ==END
TODO==

The new infrastructure has been designed with extensibility and ease-of-use in
mind.  From the perspective of an Abinit developer, adding support for the
Yaml-based approach requires two steps (??):

1. Implement the output of the YAML document in Fortran using the pre-existent
   API.  Associate a tag and possibly a label to the new document.
2. Modify the python code in pymods to associate a python object to the Yaml
   document

An example will help clarify.  Let's assume we want to implement the output of
the `!ETOT` dictionary with the different components of the total free energy.
In Fortran, we will have to write something like:

```fortran
use m_neat

call foo()
call bar()
```

To connect the new YAML document to the python infrastructure, we have 
to modify `~pymods/yaml_tools/structures.py` by adding the following lines:

```python
@yaml_auto_map  # register the following class as a known tag
                # and give it a basic set of methods so it will save
                # all fields from the YAML data as attributes
class Etot(object):
    __yaml_tag = 'ETOT'

    def __init__(self, label='nothing', comment='no comment'):
        self.label = label
        self.comment = comment

    @classmethod
    def from_map(cls, map):
        new = super(Etot, cls).from_map(map)
        new.components = {  # we want to have access to components as a coherent
                            # list
            name: value for name, value in new.__dict__.items()
            if name not in [
                'Etotal',
                'label',
                'comment',
                'Band energy',
                'Total energy(eV)'
            ]
        }
        return new
```

This code registers the object *Etot* in the Pyyaml library and instructs the Yaml parser to instantiate
this class when a document with the *!ETOT* tag is encountered.
For further details about the PyYaml API please consult the [official documentation](https://pyyaml.org/wiki/PyYAMLDocumentation).
      

## Test specification draft

### Declaration from ABINIT input file

In the standard approach, the parameters governing the execution of the test
are specified in the `TEST_INFO` section located at the end of the file.
The options are given in the [INI file format](https://en.wikipedia.org/wiki/INI_file).
The integration of the new YAML-based tests with the pre-existent infrastructure 
is obtained via two modifications of the current specifications.
More specifically:

- the *files_to_test* section now accepts the optional argument *use_yaml*. 
  Possible values are:
  
    * "yes" --> activate YAML based test (default)
    * "no" -->  do not use YAML test
    * "only" --> use YAML based test but not legacy fldiff algorithm

- a new (optional) section `[yaml_test]` is added with two possible options: 

    * *test* --> the value is used as the YAM test specification YAML
    * *file* --> path to the file containing the test specification

### Test specification syntax

To explain how to build a test specification I will go step-by-step through
the definition of the _paral[86]_ test.

First we have to introduce the main concepts.

Data tree
: The processed data is a YAML document. Therefore it can be seen as a tree
  whose nodes have a label and leaf are scalars or special data structures
  identified by a tag (not all tags mark a leaf)

Config tree
: The configuration also take the form of a tree. Its nodes are
  _specializations_ and its leaf are _parameters_ or _constraints_.
  Its structure match the structure of the _data tree_ so one can define rules
  (constraint and parameters) that apply to a specific place of the _data tree_.

Specialization
: The rules defined under a specialization will apply only on the matching node
  of the _data tree_ and its children.

Constraint
: A constraint is a condition one impose for the test to succeed.

Parameter
: A parameter is a value that can be used by constraints to modify their
  behavior.

!!! tip

    To get the list of constraints and parameters, run:
    
        ~abinit/tests/testtools.py explore
    
    and type `show *`. You can then type for example `show tol_eq` to learn more
    about a specific constraint or parameter.

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

The `tol_abs` keyword is a _constraint_. It will check for each children of
`Etot` that the value doesn't differ of more than 1.0e-7 Ha. If one does, the
test will fail and the report will tell you what component is wrong.

In the `Etot` document one can find the total energy in eV. Since the unit is
different the absolute tolerance have not the same impact on the precision.
We want to achieve the same relative precision on this term but we cannot
achieve the same absolute precision. Though we will __specialize__ the rules
and use the `tol_rel`. The name of the field holding the total energy in eV in 
the `Etot` document is `Total energy (eV)`. We could do something like this:

```yaml
Etot:
    tol_abs: 1.0e-7
    Total energy (eV):
        tol_rel: 1.0e-10
```

However the `tol_abs` constraint defined in `Etot` is _inherited_ by `Total
energy (eV)` which mean we will not only apply `tol_rel` with a tolerance of
1.0e-10 but also `tol_abs` with a tolerance of 1.0e-7. Most of the time it is
what we need even though here it's not. So we will just use a different
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
with real values. We now will have a look at the `results_gs` document. It come
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

Here we have a more examples of what can be found in a document. We have real
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

For simplicities sake I will only write the `results_gs` part in next examples.
First we want to check the integers and real values. This is what the `tol_rel`
is for.
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

    Within the explore shell `show ceil` will tell you among other informations
    what are the constraints disabled by the use of `ceil`.

Fields with the `!Tensor` tags are leaf of the tree. Then the tester routine
won't try to compare each individual coefficient with `tol_rel`. However we
still want to check that it does not change too much. For that we use the
`tol_vec` constraint which apply to all arrays derived from Numpy arrays (most
arrays with a tag). `tol_vec` check the euclidian norm of the differences
between reference and tested array. Since we also want to apply this constraint
to `cartesian_force` we will define the constraint a the top level of
`results_gs`.

```yaml
results_gs:
    tol_rel: 1.0e-8
    convergence:
        ceil: 3.0e-7
    tol_vec: 1.0e-5
```



==TODO==
Continue the How To with the explanations of the filters
==END TODO==


In some cases it is important to define different behavior depending on the
state of iterations (dtset, image...).  This is possible thanks to the so-called __filters__. 
Filters allow users to add special constraints and parameters when treating documents matching a given set of dtset, image etc. 
To define a filter, the developer uses the special node _filters_, each child of this node
is a filter.  The label of the child defines the name of the filter and its
children defines the set it matches. The below example uses two filters that
simply match one specific dataset.

```yaml
lda:
    results_gs:
        tol_abs: 1.0e-6

dmft:
    results_gs:
        cartesian forces:
            ignore

filters:
    lda:
        dtset: 1
    dmft:
        dtset: 2
```

A filter can specify all currently known iterators: dtset, timimage, image, and time. 
For each iterator a set of integers can be defined with three methods:

- a single integer value
- a YAML list of values
- a mapping with the optional members "from" and "to" specifying the boundaries (both
  included) of the integer interval. If "from" is omitted, the default is 1. If
  "to" is omitted the default is no upper boundary.

Several filters can be used for the same document if they overlap. However, it
is fundamental that an order of specificity can be determined which means that
overlapping filters must absolutely be included in each other. Though the first
example below is fine because _f2_ is included in _f1_ but the second example
will raise an error because _f4_ is not included in _f3_.

```yaml
# this is fine
filters:
    f1:
        dtset:
            from: 2
            to: 7
        image:
            from: 4

    f2:
        dtset: 7
        image:
        - 4
        - 5
        - 6
```

```yaml
# this will raise an error
filters:
    f3:
        dtset:
            from: 2
            to: 7
        image:
            from: 4

    f4:
        dtset: 7
        image:
            from: 1
            to: 5
```

When a test is defined, the default tree is overridden by the user defined
tree. When a filtered tree is used it overrides the less specific tree. As you
can see, the overriding process is used everywhere. Though, it is important to
know how it works. By default, only what is explicitly specified is overridden
which means that if a constraint is defined at a deeper level on the default
tree than what is done on the new tree, the original constraints will be kept.
For example let `f1`  and `f2` two filters such that `f2` is included in `f1`.

```yaml
f1:
    results_gs:
        tol_abs: 1.0e-6
        convergence:
            ceil: 1.0e-6
            diffor:
                1.0e-4

f2:
    results_gs:
        tol_rel: 1.0e-7
        convergence:
            ceil: 1.0e-7

filters:
    f1:
        dtset: 1
    f2:
        dtset: 1
        image: 5
```

When the tester will reach the fifth image of the first dataset, the config tree
used will be the following:

```yaml
results_gs:
    tol_abs: 1.0e-6
    tol_rel: 1.0e-7  # this has been appended without modifying anything else
    convergence:
        ceil: 1.0e-7  # this one have been overridden
        diffor:
            1.0e-4  # this one have been kept
```

If this is not the behavior you need, you can use the "hard reset marker".
Append `!` to the name of the specialization you want to override to completely
replace it. Let the `f2` tree be:

```yaml
f2:
    results_gs:
        convergence!:
            ceil: 1.0e-7
```

and now the resulting tree for the fifth image of the first dataset is:

```yaml
results_gs:
    tol_abs: 1.0e-6
    convergence:  # the whole convergence node have been overriden
        ceil: 1.0e-7
```

!!! tip

    Here again the `explore` shell could be of great help to know what is inherited
    from the other trees and what is overridden.


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

:   This tool allows the user to *explore* a test configuration file. It provides a
    shell like interface in which the user can move around the tree of the
    configuration file and print the constraints defined by the test. It also provides
    documentation about constraints and parameters via the *show* command. 

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
should help users to grasp the meaning of this file.

The third is the file *~abinit/tests/pymods/yaml_tests/structures.py*. It
defines the structures used by the YAML parser when encountering a tag (starting with !), 
or in some cases when reaching a given pattern (__undef__ for example). Even
if the abstraction layer on top of the _yaml_ module should help, it is better to
have a good understanding of more "advanced" python concepts like _inheritance_,
_decorators_, _classmethod_ etc.

### Code structure

As long as possible, it is better to include any extensions to the existing code
in these three entry points to keep the system consistent. However, this system
will eventually prove to be limited in some way. Hence, here is the general
structure of the system. (MG: Please clarify this point)

The whole system consists of three main parts:

Fldiff algorithm

:   *~abinit/tests/pymods/fldiff.py* implements the legacy algorithm and the creation of the final report.
    This module represents the interface between the yaml specific tools and the legacy test suite tools. 


Interface with Pyyaml library

:  *~abinit/tests/pymods/data_extractor.py*,
   *~abinit/tests/pymods/yaml_tools/\_\_init\_\_.py* and
   *~abinit/tests/pymods/yaml_tools/structures.py* constitute the Abinit output
   parser. *data_extractor.py* identify and extract the YAML documents from the
   source, *__init__.py* provide generic tools base on pyyaml to parse the
   documents and *structures.py* provide the classes that are used by YAML to
   handle tags. *~abinit/tests/pymods/yaml_tools/register_tag.py* define the
   abstraction layer used to simplify the code in *structures.py*. It directly
   deals with PyYaml black magic.


Yaml parsers and tools

: the other files in *~abinit/tests/pymods/yaml_tools* are dedicated to the
  testing procedure. *meta_conf_parser.py* provide the tools to read and
  interpret the YAML configuration file, namely __ConfParser__ the main parsing
  class, __ConfTree__ the tree configuration "low-level" interface (only handle
  a single tree, no access to the inherited properties, etc...) and __Constraint__
  the representation of a specific test in the tree.  *conf_parser.py* use
  __ConfParser__ to register actuals constraints and parameters. It is the place
  to define actual tests. *driver_test_conf.py* define the high level access to
  the configuration, handling several trees, applying filters etc. Finally,
  *tester.py* is the main test driver. It browses the data tree and use the
  configuration to run tests on Abinit output. It produces an __Issue__ list that
  will be used by *fldiff.py* to produce the report.

### Add a new parameter
To have the parser recognise a new token as a parameter one should edit the
*~abinit/tests/pymods/conf_parser.py*.

The `conf_parser` variable in this file have a method `parameter` to register a
new parameter. Arguments are the following:

- `token`: mandatory, the name used in configuration for this parameter
- `default`: optional (`None`), the default value used when the parameter is
  not available in the configuration.
- `value_type`: optional (`float`), the expected type of the value found in
  the configuration.
- `inherited`: optional (`True`), whether or not an explicit value in the
  configuration should be propagated to deeper levels.

### Adding a constraint
To have the parser recognise a new token as a constraint one should also edit the
*~abinit/tests/pymods/conf_parser.py*.

`conf_parser` have a method `constraint` to register a new constraint It is
supposed to be used as a decorator that takes arguments. The arguments are all
optional and are the following:

- `value_type`: (`float`) the expected type of the value found in the
  configuration.
- `inherited`: (`True`), whether or not the constraint should be propagated to
  deeper levels.
- `apply_to`: (`'number'`), the type of data this constraint can be applied to.
  This can be a type or one of the special strings (`'number'`, `'real'`,
  `'integer'`, `'complex'`, `'Array'` (refer to numpy arrays), `'this'`) `'this'`
  is used when the constraint should be applied to the structure where it is
  defined and not be inherited.
- `use_params`: (`[]`) a list of names of parameters to be passed as argument to
  the test function.
- `exclude`: (`set()`) a set of names of constraints that should not be applied
  when this one is.

The decorated function contains the actual test code. It should return `True` if
the test succeed and either `False` or an instance of `FailDetail` if the test
failed. If the test is simple enougth one should use `False`.  However if the
test is compound of several non-trivial checks `FailDetail` come in handy to
tell the user which part failed. Use `FailDetail('some explanations')` to
provide details about the failure.


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
  introduced in trunk/develop when important dimension are increased (e.g. [[chksymbreak]]  > 0 with [[ngkpt]] > 30 **3).
  Other possible applications: monitor the memory allocated to detect possible regressions.

Interface with AbiPy and Abiflows

: AbiPy has its own set of integration tests but here we mainly focus on the python layer
  without testing for numerical values. 
  Still it would be nice to check for numerical reproducibility, especially when it comes to 
  workflows that are already used in production for high-throughput applications (e.g. DFPT). 
