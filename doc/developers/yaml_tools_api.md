---
description: Tools available for YAML based tests
authors: TC
---

# YAML-based test suite: Python tools and API
## Code structure

This part aims at giving an overview of the code structure to make easier future
developments.

We can identify three logical units in the code:
- The fldiff algorithm
- The interface with PyYAML and the parsing of output data
- The parsing of test configuration and the testing logic

### Fldiff algorithm

The extraction of data from the files is performed in
*~abinit/tests/pymods/data_extractor.py*. It take care of the identification of
YAML documents in the bulk of the file and also handle the meta-character logic
of the original fldiff.pl.

*~abinit/tests/pymods/fldiff.py* is the main driver called by the *testsuite.py*
function. It implement the legacy fldiff algorithm, call the *yaml_tools*
modules to perform the YAML based tests and produce the final report.


### Interface with PyYAML library

The use of PyYAML library is abstracted away mainly by two modules of *~abinit/tests/pymods/yaml_tools/*.

In *__init__.py* are defined:
- the `yaml_parse` function that can parse documents from Abinit output files
  with all registered tag available.
- the Document class that give an interface to an extracted document with its
  metadata. The document is parsed with `yaml_parse` on demand.
- flags `is_availabe` (`True` if both PyYAML and Numpy are available) and
  `has_pandas` (`True` if Pandas is available)

In *register_tag* are defined tools to abstract the tag registration process
(see below for further details).

### YAML-based testing logic

The logic of parsing a configuration file is in *meta_conf_parser.py*. It
contains the creation of config trees, the logic to register constraints and
parameters and the logic to apply constraints.

The `Constraint` class hosts the code to build a constraint, to identify candidates
to its application and to apply it. The identify is base on the type of the
candidate which is compared to a reference type or set of types.

## Tag registration and classes implicit methods and attributes

### Basic tag registration tools

*register_tag.py* provide several decorators to be applied to classes defining
structures.

`yaml_map`
: Register a structure based on a YAML mapping. The class have to provide a
  method `from_map` that take a `dict` in argument and return an instance of
  the class. `from_map` should be a classmethod, but might be a normal method if
  the constructor accept 0 arguments.

`yaml_seq`
: Register a structure based on a YAML sequence/list. The class have to provide
  a `from_seq` method that take a `list` in argument and return an instance of the
  class.

`yaml_scalar`
: Register a structure based on a YAML scalar/anything that does not fit in the
  two other types but can be a string. The class have to provide a `from_scalar`
  method that take a `str` in argument and return an instance of the class. It
  is useful to have complex parsing. A practical case is the tables parsed by
  pandas.

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
when available. They are never required but can change some behaviour when
provided.

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
: Standart python implicit method for iterables. If an object has it the test
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

### Actual tag registration

The tag registrations 

## Constraint and parameters registration

