#!/usr/bin/env python

try:
    import yaml
except ImportError:
    raise ImportError("pyyaml package is not installed. Install it with `pip install pyyaml`")

with open('characteristics.yml', 'r') as f:
    list_chars = yaml.load(f)

with open('sections.yml', 'r') as f:
    list_sections = yaml.load(f)

list_specials = [
    ('AUTO_FROM_PSP', 'Means that the value is read from the PSP file'),
    ('CUDA', 'True if CUDA is enabled (compilation)'),
    ('ETSF_IO', 'True if ETSF_IO is enabled (compilation)'),
    ('FFTW3', 'True if FFTW3 is enabled (compilation)'),
    ('MPI_IO', 'True if MPI_IO is enabled (compilation)'),
    ('NPROC', 'Number of processors used for Abinit'),
    ('PARALLEL', 'True if the code is compiled with MPI'),
    ('SEQUENTIAL', 'True if the code is compiled without MPI'),
]


class literal(str): pass


def literal_unicode_representer(dumper, data):
    return dumper.represent_scalar('tag:yaml.org,2002:str', data, style='|')


yaml.add_representer(literal, literal_unicode_representer)


class Variable(yaml.YAMLObject):
    vartype = ''  # String containing the type
    characteristic = None  # String containing the characteristics
    definition = None  # String containing the mnemonics
    dimensions = None  # Array containing either int, formula or another variable
    defaultval = None  # Either constant number, formula or another variable
    text = None  # Description (str)
    varname = None  # Name of the variable (str)
    commentdefault = None
    commentdims = None
    section = None
    range = None
    requires = None
    excludes = None

    yaml_tag = u'!variable'

    def attrs(self):
        return ['vartype', 'characteristic', 'definition', 'dimensions', 'defaultval', 'text',
                'varname', 'section']

    def __init__(self, vartype=None, characteristic=None,
                 definition=None, dimensions=None, default=None,
                 text=None, varname=None, section=None, range=None,
                 commentdefault=None, commentdims=None):
        self.vartype = vartype
        self.characteristic = characteristic
        self.definition = definition
        self.dimensions = dimensions
        self.defaultval = default
        self.text = literal(text)
        self.varname = varname
        self.section = section
        self.commentdefault = commentdefault
        self.commentdims = commentdims
        self.range = range

    @classmethod
    def from_array(cls, array):
        return Variable(vartype=array["vartype"], characteristic=array["characteristic"],
                        definition=array["definition"], dimensions=array["dimensions"],
                        default=array["default"], text=array["text"], varname=array["varname"],
                        section=array["section"], range=array["range"], commentdefault=array["commentdefault"],
                        commentdims=array["commentdims"])

    def __str__(self):
        return "Variable " + str(self.varname) + " (default = " + str(self.defaultval) + ")"


class ValueWithUnit(yaml.YAMLObject):
    value = None
    units = None
    yaml_tag = u'!valuewithunit'

    def __init__(self, value=None, units=None):
        self.value = value
        self.units = units

    def __str__(self):
        return str(self.value) + " " + str(self.units)

    def __repr__(self):
        return str(self)


def valuewithunit_representer(dumper, data):
    return dumper.represent_mapping('!valuewithunit', data.__dict__)


class Range(yaml.YAMLObject):
    start = None
    stop = None

    yaml_tag = u'!range'

    def __init__(self, start=None, stop=None):
        self.start = start
        self.stop = stop

    def isin(self, value):
        isin = True
        if self.start is not None:
            isin = isin and (self.start <= self.value)
        if stop is not None:
            isin = isin and self.stop > self.value
        return str(self)

    def __repr__(self):
        if self.start is not None and self.stop is not None:
            return "[" + str(self.start) + " .. " + str(self.stop) + "]"
        if self.start is not None:
            return "[" + str(self.start) + "; ->"
        if self.stop is not None:
            return "<-;" + str(self.stop) + "]"
        else:
            return None


class ValueWithConditions(yaml.YAMLObject):
    yaml_tag = u'!valuewithconditions'

    def __repr__(self):
        s = ''
        for key in self.__dict__.keys():
            if key != 'defaultval':
                s += str(self.__dict__[key]) + ' if ' + str(key) + ',\n'
        s += str(self.defaultval) + ' otherwise.\n'
        return s

    def __str__(self):
        return self.__repr__()


class MultipleValue(yaml.YAMLObject):
    number = None
    value = None

    yaml_tag = u'!multiplevalue'

    def __init__(self, number=None, value=None):
        self.number = number
        self.value = value

    def __repr__(self):
        if self.number is None:
            return "*" + str(self.value)
        else:
            return str(self.number) + "*" + str(self.value)
