---
description: Tools available for YAML based tests
authors: TC
---

# YAML-based test suite: Fortran and Python API

This page gives an overview of the internal implementation in order 
to facilitate future developments and extensions.
In the first part we discuss how to write new Yaml documents using the low-level Fortran API.
In the second part, we describe the python implementation and the steps required to 
register a new Yaml document in the framework.
Finally, we provide recommendations for the name of the variables used in Yaml documents in 
order to improve readability and facilitate the integration between Yaml documents and the python code.

## Fortran API

The low-level Fortran API consists of the following modules:

**m_pair_list.F90**:

: a Fortran module for constructing dictionaries mapping strings to values.
  Internally, the dictionary is implemented in C in terms of a list of key-value pairs.
  The pair list can dynamically hold either real, integer or string values.
  It implements the basic getter and setter methods and an iterator system to allow looping
  over the key-value pairs easily in Fortran.
  <!--
  Lists have been chosen over hash tables because, in this particular context,
  performance is not critical and the O(n) cost required to locate a key is negligible
  given that there will never be a lot of elements in the table.
  -->

**m_stream_string.F90**:

:  a Fortran module providing tools to manipulate variable length
  strings with a file like interface (**stream** object)

**m_yaml_out.F90:**

:  a Fortran module to easily output YAML documents from Fortran.
  This module represents the main entry point for client code and may be used
  to implement methods to output Fortran objects.
  It provides basic routines to produce YAML documents based on a mapping at the root
  and containing well formatted (as human readable as possible and valid YAML) 1D and 2D arrays of numbers
  (integer or real), key-value mapping (using pair lists), list of key-value
  mapping and scalar fields.
  Routines support tags for specifying special data structures.

**m_neat.F90** (NEw Abinit Test system):

: a higher level module providing Fortran procedures to create specific
  documents associated to important physical properties. Currently routines
  for total energy components and ground state results are implemented.

<!--
A more detailed example is provided [below](#creating-a-new-document-in-fortran).
The new infrastructure has been designed with extensibility and ease-of-use in mind.
-->
From the perspective of the developer, adding a new YAML document requires three steps:

1. Implement the output of the YAML document in Fortran using the pre-existent API.
   Associate a **unique tag** to the new document. The tag identify the document and
   its structure.

2. register the tag and the associated class in the python code. 
   Further details about this procedure can be found below.

3. Create a YAML configuration file defining the quantities that should be compared
   with the corresponding tolerances.

An example will help clarify. Let's assume we want to implement the output of
the `!Etot` document with the different components of the total free energy.
<!--
The workflow is split in two parts. The idea is to separate the computation
from the composition of the document, and separate the composition of the
document from the actual rendering.
-->
The Fortran code that builds the list of (key, values) entries will look like:
**MG: The Fortran examples should be refactored and merged with the text in Creating a new document in Fortran**

```fortran
! Import modules
use m_neat, only: neat_etot
use m_pair_list, only: pair_list

real(dp) :: etot
type(pair_list) :: e_components  ! List of (key, value) entries

call e_components%set("comment", s="Total energie and its components")
call e_components%set("Total Energy", r=etot)

! footer of the routine, once all data have been stored
call neat_etot(e_components, unit)
```

In *47_neat/m\_neat.F90* we will implement *neat\_etot*:

```fortran
subroutine neat_etot(components, unit)

type(pair_list), intent(in) :: components
integer, intent(in) :: unit

call yaml_single_dict("Etot", "", components, 35, 500, width=20, &
                      file=unit, real_fmt='(ES20.13)')
 !
 ! 35 -> max size of the labels, needed to extract from the pairlist,
 ! 500 -> max size of the strings, needed to extract the comment
 !
 ! width -> width of the field name side, permit a nice alignment of the values

end subroutine neat_etot
```

### Creating a new document in Fortran

Developers are invited to browse the sources of `m_neat` and `m_yaml_out` to
have a comprehensive overview of the available tools. Indeed the routines are
documented and commented and creating a hand-written reference is likely to go
out of sync quicker than on-site documentation. Here we show a little example to
give the feeling of the process.

The simplest way to create a YAML document have been introduced above with the
use of `yaml_single_dict`. However this method is limited to scalars only. It is
possible to put 1D and 2D arrays, dictionaries and even tabular data.
The best way to create a new document is to have a routine `neat_my_new_document` 
that will take all required data in argument and will do all the formatting work at once.
This routine will declare a `stream_string` object to build the YAML document
inside, open the document with `yaml_open_doc`, fill it, close it with
`yaml_close_doc` and finally print it with `wrtout_stream`.
Here come a basic skeleton of `neat` routine:

```fortran
subroutine neat_my_new_document(data_1, data_2,... , iout)
  ! declare your pieces of data... (arrays, numbers, pair_list...)
  ...
  integer,intent(in) :: iout  ! this is the output file descriptor
!Local variables-------------------------------
  type(stream_string) :: stream

  ! open the document
  call yaml_open_doc('my label', 'some comments on the document', stream=stream)

  ! fill the document
  ...

  ! close and output the document
  call yaml_close_doc(stream=stream)
  call wrtout_stream(stream, iout)

end subroutine neat_my_new_document
```

Suppose we want a 2D matrix of real number with the tag `!NiceMatrix` in our document
for the field name 'that matrix' we will add the following to the middle section:

```fortran
call yaml_add_real2d('that matrix', dimension_1, dimension_2, mat_data, &
                     tag='NiceMatrix', stream=stream)
```

Other `m_yaml_out` routines provide a similar interface: first the label, then the
data and its structural metadata, then a bunch of optional arguments for
formatting, adding tags, tweaking spacing etc.

## Python API

<!--
Roughly speaking, one can group the python modules into three different parts:

- the *fldiff* algorithm
- the interface with the PyYAML library and the parsing of the output data
- the parsing of the YAML configuration file and the testing logic
-->

- The fldiff algorithm has been slightly modified to extract YAML
  documents from the output file and store them for later treatment.

- Newt tools have been written to facilitate the creation of new python classes
  corresponding to YAML tags for implementing new logic operating on the extracted data.
  These tools are very easy to use via class decorators.

- These tools have been used to create basic classes for futures tags, among
  other classes that directly convert YAML list of numbers into NumPy arrays.
  These classes may be used as examples for the creation of further tags.

- A parser for test configuration have been added and all facilities to do
  tests are in place.

- A command line tool `testtools.py` to allow doing different manual actions
  (see Test CLI)


### fldiff algorithm

The machinery used to extract data from the output file is defined in *~abinit/tests/pymods/data_extractor.py*. 
This module takes cares of the identification of the YAML documents inside the output file including 
the treatment of the meta-characters found in the first column following the original fldiff.pl implementation.

*~abinit/tests/pymods/fldiff.py* is the main driver called by *testsuite.py*.
This part implements the legacy fldiff algorithm and uses the objects defined in *yaml_tools*
to perform the YAML-based tests and produce the final report.

### Interface with the PyYAML library

The interface with the PyYAML library is implemented in *~abinit/tests/pymods/yaml_tools/*.
The `__init__.py` module exposes the public API that consists in:

- the `yaml_parse` function that parses documents from the output files
- the `Document` class that gives an interface to an extracted document with its metadata.

<!--
- flags `is_available` (`True` if both PyYAML and Numpy are available) and
  `has_pandas` (`True` if Pandas is available)
-->

- *register_tag.py* defines tools to register tags.
   See below for further details.

### YAML-based testing logic

The logic used to parse configuration files is defined in *meta_conf_parser.py*. 
It contains the creation of config trees, the logic to register constraints and
parameters and the logic to apply constraints.
The `Constraint` class hosts the code to build a constraint, to identify candidates
to its application and to apply it. The identification is based on the type of the
candidate which is compared to a reference type or as set of types.

## How to extend the test suite

### Main entry points

There are three main entry points to the system that are discussed in order of increasing complexity.
The first one is the YAML configuration file.
It does not require any Python knowledge, only a basic comprehension of the conventions
used for writing tests. Being fully *declarative* (no logic) it should be quite easy
to learn its usage from the available examples.
The second one is the *~abinit/tests/pymods/yaml_tests/conf_parser.py* file. It
contains the declarations of the available constraints and parameters. A basic
python understanding is required in order to modify this file. Comments and doc strings
should help users to grasp the meaning of this file. More details available
[in this section](./yaml_tools_api#constraints-and-parameters-registration).
The third one is the file *~abinit/tests/pymods/yaml_tests/structures/*. It
defines the structures used by the YAML parser when encountering a tag (starting
with !), or in some cases when reaching a given pattern (__undef__ for example).
The *structures* directory is a package organised by features (ex: there is a
file *structures/ground_state.py*). Each file define structures for a given
feature. All files are imported by the main script `structures/__init__.py`.
Even if the abstraction layer on top of the _yaml_ module should help, it is
better to have a good understanding of more "advanced" python concepts like
*inheritance*, *decorators*, *classmethod* etc.

## Tag registration and classes implicit methods and attributes

### Basic tag registration tools

The *register_tag.py* module defines python decorators that are used to build python classes 
associated to Yaml documents.

`yaml_map`
: Registers a structure based on a YAML mapping. The class must provide a
  method `from_map` that receives a `dict` object as argument and returns an instance of
  the class. `from_map` should be a class method, but might be a normal method if
  the constructor accepts zero arguments.

`yaml_seq`
: Register a structure based on a YAML sequence/list. The class must provide
  a `from_seq` method that takes a `list` as input and returns an instance of the class.

`yaml_scalar`
: Register a structure based on a YAML scalar/anything that does not fit in the
  two other types but can be a string. The class must provide a `from_scalar`
  method that takes a `string` in argument and returns an instance of the class. It
  is useful to have complex parsing. A practical case is the CSV table parsed by pandas.

`yaml_implicit_scalar`
: Provide the possibility to have special parsing without explicit tags. The
  class it is applied to should meet the same requirements than for `yaml_scalar`
  but should also provide a class attribute named `yaml_pattern`. It can be either
  a string or a compile regex and it should match the expected structure of the
  scalar.

`auto_map`
: Provide a basic general purpose interface for map objects:
  - dict-like interface inherited from `BaseDictWrapper` (`get`, `__contains__`,
    `__getitem__`, `__setitem__`, `__delitem__`, `__iter__`, `keys` and `items`)
  - a `from_map` method that accept any key from the dictionary and add them as
  attributes.
  - a `__repr__` method that show all attributes

`yaml_auto_map`
: equivalent to `auto_map` followed by `yaml_map`

`yaml_not_available_tag`
: This is a normal function (not a decorator) that take a tag name and a message as
  arguments and will register the tag but will raise a warning each time the tag
  is used. An optional third argument `fatal` replace the warning by an error if
  set to `True`.


### Implicit methods and attributes

Some class attributes and methods are used here and there in the tester logic
when available. They are never required but can change some behaviour when provided.

`_is_base_array` Boolean class attribute
: Assert that the object derived from `BaseArray` when set to `True` (`isinstance`
  is not reliable because of the manipulation of the `sys.path`)

`_not_available` Boolean class attribute
: set to `True` in object returned when a tag is registered with
  `yaml_not_available_tag` to make all constraints fail with a more
  useful message.

`is_dict_like` Boolean class attribute
: Assert that the object provide a full `dict` like interface when set to
  `True`. Used in `BaseDictWrapper`.

`__iter__` method
: Standard python implicit method for iterables. If an object has it the test
  driver will crawl its elements.

`get_children` method
: Easy way to allow the test driver to crawl the children of the object.
  It does not take any argument and return a dictionary `{child_name: child_object ...}`.

`has_no_child`  Boolean class attribute
: Prevent the test driver to crawl the children even if it has an `__iter__`
  method (not needed for strings).

`short_str` method
: Used (if available) when the string representation of the object is too long
  for error messages. It should return a string representing the object and not
  too long. There is no constraints on the output length but one should keep
  things readable.

## Constraints and parameters registration

### Add a new parameter

To have the parser recognise a new token as a parameter, one should edit the
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

Example:

```python
conf_parser.parameter('tol_eq', default=1e-8, inherited=True)
```

### Adding a constraint

To have the parser recognise a new token as a constraint, one should also edit
*~abinit/tests/pymods/conf_parser.py*.

`conf_parser` has a method `constraint` to register a new constraint. It is
supposed to be used as a decorator (on a function) that takes keywords arguments. 
The arguments are all optional and are the following:

- `name`: (`str`) the name to be used in config files. If not
  specified, the name of the function is used.
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
- `handle_undef`: (`True`) whether or not the special value `undef`
  should be handled before calling the test function.
  If `True` and a `undef` value is present in the data the test will fail or
succeed depending on the value of the special parameter `allow_undef`.
  If `False`, `undef` values won't be checked. They are equivalent to *NaN*.

The decorated function contains the actual test code. It should return `True` if
the test succeed and either `False` or an instance of `FailDetail` (from
*common.py*) if the test failed.

If the test is simple enough one should use `False`.  However if the
test is compound of several non-trivial checks `FailDetail` come in handy to
tell the user which part failed. When you want to signal a failed test and
explaining what happened return `FailDetail('some explanations')`. The message
passed to `FailDetail` will be transmitted to the final report for the user.

Example with `FailDetail`:

```python
@conf_parser.constraint(exclude={'ceil', 'tol_abs', 'tol_rel', 'ignore'})
def tol(tolv, ref, tested):
    '''
        Valid if both relative and absolute differences between the values
        are below the given tolerance.
    '''
    if abs(ref) + abs(tested) == 0.0:
        return True
    elif abs(ref - tested) / (abs(ref) + abs(tested)) >= tolv:
        return FailDetail('Relative error above tolerance.')
    elif abs(ref - tested) >= tolv:
        return FailDetail('Absolute error above tolerance.')
    else:
        return True
```

### How to add a new tag

Pyyaml offer the possibility to directly convert YAML documents to a
Python class using tags. To register a new tag, edit the file `~abinit/tests/pymods/yaml_tools/structures.py`. 
In this file there are several
classes that are *decorated* with one of `@yaml_map`, `@yaml_scalar`, `@yaml_seq` or `@yaml_auto_map`.
These [decorators](https://realpython.com/primer-on-python-decorators/) are the functions that actually
register the class as a known tag.
Whether you should use one or another depend on the organization of the data in the YAML document.
Is it a dictionary, a scalar or a list/sequence ?

In the majority of the cases, one wants the `tester` to browse the children of the structure and apply
the relevant constraints on it and to access attributes through their original name
from the data. In this case, the simpler way to register a new tag is to use
`yaml_auto_map`. For example to register a tag Etot that simply register all
fields from the data tree and let the tester check them with `tol_abs` or
`tol_rel` we would put the following in *~abinit/tests/pymods/yaml_tools/structures/ground_state.py*:

```python
@yaml_auto_map
class Etot(object):
    pass
```

Now the name of the class "Etot" is associated with a tag recognized by the YAML
parser.

`yaml_auto_map` does several things for us:

- it gives the class a `dict` like interface by defining relevant methods, which
  allows the `tester` to browse the children
- it registers the tag in the YAML parser
- it automatically registers all attributes found in the data tree as attributes
  of the class instance. These attributes are accessible through the attribute
  syntax (ex: `my_object.my_attribute_unit`) with a normalized name (basically remove
  characters that cannot be in a python identifier like spaces and punctuation characters)
  or through the dictionary syntax with their original name (ex: `my_object['my attribute (unit)']`)

Sometimes one wants more control over the building of the class instance.
This is what `yaml_map` is for.
Let suppose we still want to register tag but we want to select only a subset of
the components for example. We will use `yaml_map` to gives us control over the
building of the instance. This is done by implementing the `from_map` class
method (a class method is a method that is called from the class instead of
being called from an instance). This method take in argument a dictionary built
from the data tree by the YAML parser and should return an instance of our class.

```python
@yaml_map
class Etot(object):
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
    # same code as the previous example

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

```yaml
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

    # same code as the previous example

    def get_children(self):
        return {
            'kin': self.kinetic,
            'hart': self.hartree,
            'xc': self.xc,
            'ew': self.ewald
        }
```

Now `tester` will be able to apply `tol_abs` and friends to the components we gave him.

If the class has a __complete__ dict-like read interface (`__iter__` yielding
keys, `__contains__`, `__getitem__`, `keys` and `items`) then it can have a
class attribute `is_dict_like` set to `True` and it will be treated as any other
node (it not longer need `get_children`). `yaml_auto_map` registered classes
automatically address these requirements.

`yaml_seq` is analogous to `yaml_map` however `to_map` became `to_seq` and the
YAML source data have to match the YAML sequence structure that is either `[a, b, c]` or

```yaml
- a
- b
- c
```

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

```python
@yaml_implicit_scalar
class YAMLComplex(complex):

    # A (quite rigid) regular expression matching complex numbers
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

## Recommend conventions for Yaml documents

This section discusses the basic conventions that should be followed when writing YAML documents in Fortran.
Note that these conventions are motivated by technical aspects that will facilitate the integration with the
python language as well as the implementation of post-processing tools.

* The tag name should use CamelCase so that one can directly map tag names to python classes
  that are usually given following this convention (see also [PEP8](https://www.python.org/dev/peps/pep-0008/)).

* The tag name should be self-explanatory.

* Whenever possible, the keywords should be a **valid python identifier**. This means one should avoid white spaces
  as much possible and avoid names starting with non-alphabetic characters. White spaces should be replaced by
  underscores. Prefer lower-case names whenever possible.

* By default, quantities are supposed to be given in **atomic units**. If the document contains quantities given
  in other units, we recommended to encode this information in the name using the syntax: `foo_eV`.

*  Avoid names starting with an underscore because these names are reserved for future additions 
   in the python infrastructure.

* `tol_abs`, `tol_rel`, `tol_vec`, `tol_eq`, `ceil`, `ignore`, `equation`, `equations`, `callback` and `callbacks` are reserved 
   keywords that shall not be used in Yaml documents.

* The *comment* field is optional but it is recommended especially when the purpose of the document is not obvious.

An example of well-formed document

```yaml
--- !EnergyTerms
comment             : Components of total free energy (in Hartree)
kinetic_energy      :  5.279019930263079807E+00
hartree_energy      :  8.846303409910728499E-01
xc_energy           : -4.035286123400158687E+00
total_energy_eV     : -2.756386620520307815E+02
...
```

An example of Yaml document that does not follow our guidelines:

```yaml
--- !ETOT
Kinetic energy      :  5.279019930263079807E+00
Hartree energy      :  8.846303409910728499E-01
_XC energy          : -4.035286123400158687E+00
Total energy(eV)    : -2.756386620520307815E+02
equation            : 1 + 1 = 3
...
```

<!--
In a nutshell:

* no savage usage of write statements!
* In the majority of the cases, YAML documents should be opened and closed in the same routine!
* To be discussed: `crystal%yaml_write(stream, indent=4)`
* `tol_rel`, `tol_abs` etc are now reserved keywords that cannot be used in YAML documents. This will confuse the parser.
  We should add a check in the Fortran code and stop the code in a reserved keyword is used. 
* One should also promote the use of python-friendly keys: white space replaced by underscore, no number as first character,
  tags should be CamelCase
* Is Yaml mode supported only for main output files? What happens to the other files in the files\_to\_test section?
<!--
Tags are now mandatory for documents
* A document with a tag is considered a standardized document i.e. a document for which there's
  an official commitment from the ABINIT community to maintain backward compatibility.
  Official YAML documents may be used by third-party software to implement post-processing tools.
-->


: The *label* field appears in all data document. It should be a unique identifier
  at the scale of an _iteration state_. The tag is not necessarily unique, it
  describes the structure of the document and there is no need to use it unless
  special logic have to be implemented for the document. 
-->

<!--
## Further developments

This new YAML-based infrastructure can be used as building block to implement the
high-level logic required by more advanced integration tests such as:

Parametrized tests

: Tests in which multiple parameters are changed either at the level of the input variables
  or at the MPI/OpenMP level.
  Typical example: running calculations with [[useylm]] in [0, 1] or [[paral_kgb]] = 1 runs
  with multiple configurations of [[npfft]], [[npband]], [[npkpt]].

Benchmarks

: Tests to monitor the scalability of the code and make sure that serious bottlenecks are not
  introduced in trunk/develop when important dimension are increased (e.g.
  [[chksymbreak]] > 0 with [[ngkpt]] > 30\*\*3).
  Other possible applications: monitor the memory allocated to detect possible regressions.

Interface with AbiPy and Abiflows

: AbiPy has its own set of integration tests but here we mainly focus on the python layer
  without testing for numerical values.
  Still it would be nice to check for numerical reproducibility, especially when it comes to
  workflows that are already used in production for high-throughput applications (e.g. DFPT).
-->
