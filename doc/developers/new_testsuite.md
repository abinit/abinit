---
description: Design of the new YAML-based test suite.
authors: TC, MG
---

# YAML-based tests: design principles and implementation details

This document discusses the new YAML-based format used to represent physical results
in the main output file.
The connection with the ABINIT test suite and the yaml syntax used to define
tolerances, parameters and constraints is also discussed.

The new infrastructure consists of a set of Fortran modules to output 
structured data in YAML format and Python code to parse the output files and analyze data.

### Motivations

## Motivations

In ABINITv8 and previous versions, the ABINIT test suite is based on input files with the
associated reference output files.
Roughly speaking, an automatic test consists in comparing the reference file with
the output file in a **line-oriented fashion**
by computing differences between floating-point numbers **without any knowledge** about
the meaning and the importance of the numerical values.
This approach is very rigid and rather limited because it obliges developers to use a
single (usually large) tolerance to account for possibly large fluctuations in
the intermediate results whereas the tolerance criteria should be ideally
applied only to the final results that are (hopefully) independent of details
such as hardware, optimization level and parallelism.

This limitation is clearly seen when comparing the results of iterative algorithms.
The number of iterations required to converge, indeed, may depend
on several factors especially when the input parameters are far from convergence
or when different parallelization schemes or stochastic methods are used.
From the point of view of code validation, what really matters is the **final converged value**
plus the **time to solution** if performance ends up being of concern.
Unfortunately, any line-by-line comparison algorithm will miserably fail in such
conditions because it will continue to insist on having the same number of
iterations (lines) in the two calculations to consider the test succeeded.
It is therefore clear that the approach used so far to validate new developments in
Abinit is not able to cope with the challenges posed by high-performance
computing and that smarter and more flexible approaches are needed to address these limitations.
Ideally, we would like to be able to

* use **different tolerances** for particular quantities that are selected by keyword
* be able to replace the check on the absolute and relative difference with a **threshold check**
  (is this quantity smaller that the give threshold?)
* have some sort of syntax to apply different rules depending on the **iteration state** e.g. the dataset index.
* execute python code (**callbacks**) that operates on the data to perform more advanced tests requiring 
  some sort of post-processing
* provide an easy-to-use **declarative interface** that allows developers to define the logic
  to compare selected quantities.

In what follows, we present this new infrastructure, the design principles and the steps required
to define YAML-based tests.

<!--
The goal of this project is to provide ABINIT developers with reusable tools to implement physics-aware tests.
With "physics-aware" we intend tests in which the developer can customize the tolerances and the logic used to
compare numerical results thanks to the fact that the new parser
is aware of the context and of the meaning of the numerical values extracted from the output file.
The new infrastructure consists of a set of Fortran modules to output
structured data in YAML format as well as Python code to parse the output files and analyze data.
One of the goals of this project is to implement a 

The new python-based infrastructure allows developers to implement more rigorous tests on the portability of
important physical quantities while ignoring intermediate results.
This goal is achieved by providing a *declarative* interface that allows developers to define
the logic used to compare selected physical quantities.
We also provide a declarative API to check that the results computed by the ab-initio code
satisfies fundamental physical and mathematical rules such as energy
conservation, Newton's third law, symmetry properties, etc.

That is not quite as simple as it sounds, especially because one of our goals is
to minimize the amount of (coding) work required to express such logic. 
There are therefore basic rules and design principles Abinit developers should be
aware of in order to take fully advantage of the new infrastructure. 
In what follows, we briefly describe the philosophy employed to implement the YAML-based
test suite, the syntax used to define tests and the modifications required to extend the framework.
-->

!!! important

    The new YAML-based testsuite relies on libraries that are not provided
    by the python standard library.
    To install these dependencies in *user* mode, use:

        pip install numpy pyyaml pandas --user

    If these dependencies are not available, the new test system will be
    disabled and a warning message is printed to the terminal.


## Implementation details

For reasons that will be clear later, implementing smart algorithms requires metadata and context.
In other words the python code needs to have some basic understanding of the meaning of the numerical values 
extracted from the output file and must be able to locate a particular property by name or by its "position"
inside the output file.
For this reason, the most important physical results are now written in the main output file (*ab_out*)
using machine-readable [YAML](https://en.wikipedia.org/wiki/YAML) documents.

A YAML document starts with three hyphens (---) followed by an
*optional* tag beginning with an exclamation mark (e.g. `!ETOT`). 
Three periods (...) signals the end of the document. 
Following these rules, one can easily write a dictionary containing the different contributions
to the total free energy using:


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

Further details about the meaning of **tags**, **labels** and their connection with the
testing infrastructure will be given in the below sections. For the time being,
it is sufficient to say that we opted for YAML because it is a *human-readable*
data serialization language already used in the log file to record important
events such as WARNINGs, ERRORs and COMMENTs (this is indeed the *protocol*
employed by AbiPy to monitor the status of Abinit calculations).
Many programming languages, including python, provide support for YAML hence
it is relatively easy to implement post processing tools
based on well-established python libraries for scientific computing such as
NumPy, SciPy, Matplotlib and Pandas. Last but not least, writing YAML in
Fortran does not represent an insurmountable problem provided one keeps the
complexity of the YAML document at a *reasonable level*.
Our Fortran implementation, indeed, supports only a subset of the YAML specifications:

* scalars
* arrays with one or two dimensions
* dictionaries mapping strings to scalars
* tables in CSV format (this is an extension of the standard)

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

<!--
In this tutorial we mainly focus on the user interface i.e. on the syntax used to define Yaml-based tests.
For the technical details related to the internal implementation see
-->

## YAML configuration file

### How to activate the YAML mode

The parameters governing the execution of the test
are specified in the `TEST_INFO` section located at the end of the input file.
The options are given in the [INI file format](https://en.wikipedia.org/wiki/INI_file).
The integration of the new YAML-based tests with the pre-existent infrastructure
is obtained via two modifications of the current specifications.
More specifically:

- the **files_to_test** section now accepts the optional argument **use_yaml**.
  The allowed values are:

    * "yes" --> activate YAML mode
    * "no" -->  do not use YAML mode (default)
    * "only" --> use YAML mode, deactivate legacy fldiff algorithm

- a new *optional* section **[yaml_test]** has been added. 
  This section contains two mutually exclusive fields:

    * *file* --> path of the YAML configuration file.
      The path is relative to the input file. A natural choice
      would be to use the same prefix as the input file e.g. "./t21.yaml" 
      is the configuration file associated to the input file "t21.in".

    * *test* --> multi-line string with the YAML specifications. This option may
      be used for really short configurations heavily relying on the default values.


An example of `TEST_INFO` section that activates the YAML mode can be found in [[test:paral_86]]:

```sh
#%%<BEGIN TEST_INFO>
#%% [setup]
#%% executable = abinit
#%% [files]
#%% psp_files = 23v.paw, 38sr.paw, 8o.paw
#%% [paral_info]
#%% nprocs_to_test = 4
#%% max_nprocs = 4
#%% [NCPU_4]
#%% files_to_test =
#%%   t86_MPI4.out, use_yaml = yes, tolnlines = 4, tolabs = 2.0e-2, tolrel = 1.0, fld_options = -easy;
#%% [extra_info]
#%% authors = B. Amadon, T. Cavignac
#%% keywords = DMFT, FAILS_IFMPI
#%% description = DFT+DMFT for SrVO3 using Hubard I code with KGB parallelism
#%% topics = DMFT, parallelism
#%% [yaml_test]
#%% file = ./t86.yaml
#%%<END TEST_INFO>
```

with the associated YAML configuration file given by:

{% dialog tests/paral/Input/t86.yaml %}

!!! important

    Until the basic document list is considered stable enough the printing of
    YAML documents is disabled by default. To enable it, add [[use_yaml]] 1 in the input file 

### Out first example of YAML configuration file 

Let's start with a minimalistic example in which we compare the components of
the total free energy in the `Etot` document with an absolute tolerance of 1.0e-7 Ha.
The YAML configuration file will look like:

```yaml
Etot:
    tol_abs: 1.0e-7
```

The `tol_abs` keyword defines the **constraint** that will applied to **all the children** of the `Etot` document.
In other words, all the entries in the `Etot` dictionary will be compared with 
an absolute tolerance of 1.0e-7 and a default value for the relative difference `tol_rel`.
In order to pass the test, all the absolute differences must be smaller that `tol_abs`.
<!--
If this condition is not fulfilled, the test will fail and python will list the entries that did not pass the test.
-->

There are however cases in which we would like to specify different tolerances for particular entries
instead of relying on the *global* tolerances.
The `Etot` document, for example, contains the total energy in eV in the `Total energy (eV)` entry.
To use a different absolute tolerance for this property, we **specialize** the rule with the syntax:

<!--
Since the unit is different, the absolute tolerance does not have the same impact on the precision.
We want to achieve the same relative precision on this term but we cannot
achieve the same absolute precision. Though we will **specialize** the rules
and use the `tol_rel` constraint. 
We could do something like this:

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
-->

```yaml
Etot:
    tol_abs: 1.0e-7
    Total energy (eV):
        tol_abs: 1.0e-5
```

<!--
Now we achieve the same relative precision and the test does not fail because
of the looser absolute precision of the total energy in eV.
-->

### Basic concepts

In the previous section, we presented a minimal example of configuration file.
In the next paragraphs we will discuss in more detail how to implement more advanced test
but before proceeding with the examples, we need to introduce some basic terminology
to facilitate the discussion.
<!--
To explain how to build a test specification we will go step-by-step through
the definition of [[test:paral_86]].
First we have to introduce the most important concepts.
-->

Document tree
: The YAML document is a dictionary that can be treated as a tree
  whose nodes have a label and leaf are scalars or special data structures
  identified by a tag (note however that not all tags mark a leaf).
  The top-level node is the YAML document.

Config tree
: The YAML configuration also takes the form of a tree where nodes are
  **specializations** and its leaf represent **parameters** or **constraints**.
  Its structure match the structure of the **document tree** thus one can define rules
  (constraint and parameters) that will be applied to a specific part of the **document tree**.

Specialization
: The rules defined under a specialization will apply only on the matching node
  of the *document tree* and its children.

Constraint
: A constraint is a condition one imposes for the test to succeed. Constraints
  can apply to leafs of the document tree or to nodes depending of the nature of the constraint.

Parameter
: A parameter is a value that can be used by the constraints to modify their behavior.

Iteration state
: An iteration state describes how many iterations of each possible level are present
  in the run (e.g. idtset = 2, itimimage = not used, image = 5, time = not used).
  It gives information on the current state of the run. Documents are implicitly
  associated to their iteration state. This information is made available to
  the test engine through specialized YAML documents using the `IterStart` tag.

!!! tip

    To get the list of constraints and parameters, run:

        ~abinit/tests/testtools.py explore

    and type `show *`. You can then type for example `show tol_eq` to learn more
    about a specific constraint or parameter.

A few conventions on documents writing

: The *label* field appears in all data document. It should be a unique identifier
  at the scale of an _iteration state_. The tag is not necessarily unique, it
  describes the structure of the document and there is no need to use it unless
  special logic have to be implemented for the document. The *comment* field is
  optional but it is recommended especially when the purpose of the document is not obvious.

### A more complicated example

The `Etot` document is the simplest possible document. It only contains fields
with real values. Now we will have a look at the `results_gs` document
that represents the results stored in the corresponding Fortran datatype used in Abinit.
The YAML document is now given by:

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

This YAML document is more complicated as it contains scalar fields, dictionaries and even 2D arrays.
**MG: Are integer values always compared without tolerance?**
Still, the parsers will be able to locate the entire document via its tag/label and address all the entries by name.
To specify the tolerance for the relative difference for all the scalar quantities in `results_gs`,
we just add a new entry to the YAML configuration file similarly to what we did for `Etot`: 

```yaml
Etot:
    tol_abs: 1.0e-7
    Total energy (eV):
        tol_abs: 1.0e-5
        tol_rel: 1.0e-10

results_gs:
    tol_rel: 1.0e-8
```

<!--
For simplicity sake I will only write the `results_gs` part in next examples.
-->
Unfortunately, such a strict value for `tol_rel` will become very problematic
when we have to compare the residues stored in the `convergence` dictionary!
In this case, it makes more sense to check that all the residues are below a certain threshold.
This is what the **ceil** constraint is for:

```yaml
results_gs:
    tol_rel: 1.0e-8
    convergence:
        ceil: 3.0e-7
```

Now the test will fail if one of the components of the `convergence` dictionary is above 3.0e-7.
Note that the `ceil` constraint automatically disables the check for `tol_rel` and `tol_abs` inside `convergence`.
In other words, all the scalar entries in `results_gs` will be compared with our `tol_rel` and the default `tol_abs`
whereas the entries in the `convergence` dictionary will be tested against `ceil`.

!!! tip

    Within the explore shell `show ceil` will list
    the constraints that are disabled by the use of `ceil` in the _exclude_ field.

Up to now we have been focusing on scalar quantities for which the concept of 
relative and absolute difference is unambiguously defined but how do we compare vectors and matrices?
Fields with the `!Tensor` tags are leafs of the tree. The tester routine
won't try to compare each individual coefficient with `tol_rel`. However we
still want to check that it does not change too much. For that purpose we use the
`tol_vec` constraint which apply to all arrays derived from `BaseArray` (most
arrays with a tag). `BaseArray` let us use the capabilities of Numpy arrays with
YAML defined arrays. `tol_vec` check the euclidean distance between the reference
and the output arrays. Since we also want to apply this constraint to
`cartesian_force`, we will define the constraint at the top level of `results_gs`.

```yaml
results_gs:
    tol_rel: 1.0e-8
    tol_vec: 1.0e-5
    convergence:
        ceil: 3.0e-7
```

### How to use filters to select documents by iteration state

Thanks to the syntax presented in the previous sections, one can customize tolerances for different 
documents and different entries.
Note however that these rules will be applied to all the documents found in the output file.
This means that we are implicitly assuming that all the different steps of the calculation
have similar numerical stability.
There are however cases in which the results of particular datasets are less numerically stable than the others.
An example will help clarify. 

The test [[test:paral_86]] uses two datasets to perform two different computations.
The first dataset computes the DFT density with LDA while
the second dataset uses the LDA density to perform a DMFT computation.
The entire calculation is supposed to take less than ~3-5 minutes hence the input parameters 
are severely under converged and the numerical noise propagates quickly through the different steps.
As a consequence, one cannot expect the DFMT results to have the same numerical stability as the LDA part. 
Fortunately, one can use **filters** to specify different convergence criteria for the two datasets.

A filter is a mechanism that allows one to associate a specific configuration to a set of **iteration states**.
A filter is defined in a separated section of the configuration file under the node `filters`.
Let's declare two filters with the syntax:

```yaml
filters:
    lda:
        dtset: 1
    dmft:
        dtset: 2
```

Here we are simply saying that we want to associate the label `lda` to all
documents created in the first dataset and the label `dmft` to all document created in the second dataset.
This is the simplest filter declaration possible.
See [here](#filters-api) for more info on filter declarations.
Now we can use our filters. First of all we will associate the configuration
we already wrote to the `lda` filter so we can have a different configuration
for the second dataset. The YAML file now reads

```yaml
filters:
    lda:
        dtset: 1
    dmft:
        dtset: 2

lda:
    Etot:
        tol_abs: 1.0e-7
        Total energy (eV):
            tol_abs: 1.0e-5
            tol_rel: 1.0e-10

    results_gs:
        tol_rel: 1.0e-8
        tol_vec: 1.0e-5
        convergence:
            ceil: 3.0e-7

```

By putting the configuration under the `lda` node, we specify that these rules
apply only to the first dataset. We will then create a new `dmft` node and create a
configuration following the same procedure as before. 
We end up with something like this:

```yaml
filters:
    lda:
        dtset: 1
    dmft:
        dtset: 2

lda:
    Etot:
        tol_abs: 1.0e-7
        Total energy (eV):
            tol_abs: 1.0e-5
            tol_rel: 1.0e-10

    results_gs:
        tol_rel: 1.0e-8
        tol_vec: 1.0e-5
        convergence:
            ceil: 3.0e-7

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
```

## Filters API

Filters provide a practical way to specify different configuration for different
states of iterations without having to rewrite everything from scratch.

### Filter declaration

A filter can specify all currently known iterators: **dtset**, **timimage**, **image**, and **time**.
For each iterator, a set of integers can be defined with three different methods:

- a single integer value (`dtset: 1`)
- a YAML list of values (`dtset: [1, 2, 5]`)
- a mapping with the optional members "from" and "to" specifying the boundaries (both
  included) of the integer interval (`dtset: {from: 1, to: 5}`). If "from" is omitted, the default is 1. If
  "to" is omitted the default is no upper boundary. 

### Filter overlapping

Several filters can apply to the same document if they overlap. However, they
are required to have a trivial order of *specificity*. Though the first example
below is fine because _f2_ is included (i.e. is more specific) in _f1_ but the
second example will raise an error because _f4_ is not included in _f3_.

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

When a test is defined, the default tree is overridden by the user-defined tree.
When a filtered tree is used, it overrides the less specific tree. Trees are
sequentially applied to the tree from the most general to the most specific one.
The overriding process is often used, though it is important to know how it
works. By default, only what is explicitly specified in the file is overridden which means
that if a constraint is defined at a deeper level on the default tree than what
is done on the new tree, the original constraints will be kept. For example let
`f1`  and `f2` be two filters such that `f2` is included in `f1`.

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
    tol_abs: 1.0e-6  # this come from application of f1
    tol_rel: 1.0e-7  # this has been appended without modifying anything else when appling f2
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

### How to use equation and callback

**equation** and **callback** are special constraints because their
actual effects are defined directly in the configuration file. They have been introduced to
increase the flexibility of the configuration file without having to change the python code.

**equation** takes a string in input. This string will
be interpreted as a python expression that must return in a number. The
absolute value of this number will be compared to the value of the `tol_eq`
parameter and if `tol_eq` is greater the test will succeed. The expression can
also result in a numpy array. In this case, the returned value if the euclidean norm of the
array that will be compared to `tol_eq` value.
A minimal example:

```yaml
Etot:
    tol_eq: 1.0e-6
    equation: 'this["Etotal"] - this["Total energy(eV)"]/27.2114'
```

**equations** works exactly the same but has a list of string as value. 
Each string is a different expression
that will be tested independently from the others. In both case the tested
object can be referred as `this` and the reference object can be referred as `ref`.


**callback** requires a bit of python coding since it will invoke a method of the
structure. Suppose we have a tag `!AtomSpeeds` associated
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

## Command line interface

The `~abinit/tests/testtools.py` script provides a command line interface to facilitate 
the creation of new tests and the exploration of the YAML configuration file.
The syntax is:

    ./testtools.py COMMAND [options]

Run the script without arguments to get the list of possible commands and use:

    ./testtools.py COMMAND --help

to display the options supported by `COMMAND`.
The list of available commands is:

fldiff

:   Interface to the *fldiff.py* module.
    This command can be used to compare output and reference files without executing ABINIT.
    It is also possible to specify the YAML configuration file with the `--yaml-conf` option so 
    that one can employ the same parameters as those used by *runtests.py*

explore

:   This command allows the user to *explore* and *validate* a YAML configuration file. 
    It provides a shell like interface in which the user can explore
    the tree defined by the configuration file and print the constraints. 
    It also provides documentation about constraints and parameters via the *show* command.
