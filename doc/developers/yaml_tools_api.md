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

## Constraints and parameters registration

### Adding a new parameter
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

Example:
```python
conf_parser.parameter('tol_eq', default=1e-8, inherited=True)
```

### Adding a constraint
To have the parser recognise a new token as a constraint one should also edit the
*~abinit/tests/pymods/conf_parser.py*.

`conf_parser` have a method `constraint` to register a new constraint It is
supposed to be used as a decorator (on a function) that takes keywords arguments. The arguments are all
optional and are the following:

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

## Filters API

Filters provide a practical way to specify different configuration for different
states of iterations without rewriting everything from scratch.

### Filter declaration

A filter can specify all currently known iterators: dtset, timimage, image, and time. 
For each iterator a set of integers can be defined with three methods:

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

When a test is defined, the default tree is overridden by the user defined tree.
When a filtered tree is used it overrides the less specific tree. Trees are
sequentially applied to the tree from the most general to the most specific.
The overriding process is often used, though it is important to know how it
works: By default, only what is explicitly specified is overridden which means
that if a constraint is defined at a deeper level on the default tree than what
is done on the new tree, the original constraints will be kept.  For example let
`f1`  and `f2` two filters such that `f2` is included in `f1`.

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
