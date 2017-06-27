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

list_topics_class = [
    ('compulsory', 'Compulsory input variables:'),
    ('basic', 'Basic input variables:'),
    ('useful', 'Useful input variables:'),
    ('internal', 'Relevant internal variables:'),
    ('prpot', 'Printing input variables for potentials:'),
    ('prfermi', 'Printing input variables for fermi level or surfaces:'),
    ('prden', 'Printing input variables for density, eigenenergies k-points and wavefunctions:'),
    ('prgeo', 'Printing input variables for geometry:'),
    ('prdos', 'Printing DOS-related input variables:'),
    ('prgs', 'Printing other ground-state input variables:'),
    ('prngs', 'Printing non-ground-state input variables:'),
    ('alch', 'Input variables related to alchemical mixing:'),
    ('job', 'Input variables for job time limits:'),
    ('slab', 'Input variables to insert a slab:'),
    ('lotf', 'Input variables for Learn On The Fly calculations:'),
    ('expert', 'Input variables for experts:'),
]

class literal(str): pass


def literal_unicode_representer(dumper, data):
    return dumper.represent_scalar('tag:yaml.org,2002:str', data, style='|')


yaml.add_representer(literal, literal_unicode_representer)


class Variable(yaml.YAMLObject):
    vartype = ''  # String containing the type
    characteristic = None  # String containing the characteristics
    mnemonics = None  # String containing the mnemonics
    dimensions = None  # Array containing either int, formula or another variable
    defaultval = None  # Either constant number, formula or another variable
    text = None  # Description (str)
    abivarname = None  # Name of the variable (str)
    commentdefault = None
    commentdims = None
    section = None
    range = None
    requires = None
    excludes = None
    topic_name = None
    topic_class = None

    yaml_tag = u'!variable'

    def attrs(self):
        return ['vartype', 'characteristic', 'mnemonics', 'dimensions', 'defaultval', 'text',
                'abivarname', 'section', 'topic_name', 'topic_class']

    def __init__(self, vartype=None, characteristic=None,
                 mnemonics=None, dimensions=None, default=None,
                 text=None, abivarname=None, section=None, range=None,
                 commentdefault=None, commentdims=None, topic_name=None, topic_class=None):
        self.vartype = vartype
        self.characteristic = characteristic
        self.mnemonics = mnemonics
        self.dimensions = dimensions
        self.defaultval = default
        self.text = literal(text)
        self.abivarname = abivarname
        self.section = section
        self.commentdefault = commentdefault
        self.commentdims = commentdims
        self.range = range
        self.topic_name = topic_name
        self.topic_class = topic_class

    @classmethod
    def from_array(cls, array):
        return Variable(vartype=array["vartype"], characteristic=array["characteristic"],
                        mnemonics=array["mnemonics"], dimensions=array["dimensions"],
                        default=array["default"], text=array["text"], abivarname=array["abivarname"],
                        section=array["section"], range=array["range"], commentdefault=array["commentdefault"],
                        commentdims=array["commentdims"], topic_name=array["topic_name"], topic_class=array["topic_class"])

    def __str__(self):
        return "Variable " + str(self.abivarname) + " (default = " + str(self.defaultval) + ")"
<<<<<<< HEAD
=======

class Topic(yaml.YAMLObject):
    topic_name = None  # String containing the "How to ?" topic name
    howto = ''     # String containing the description of the topics, to be echoed after "How to" ...

    yaml_tag = u'!topic'

    def attrs(self):
        return ['topic_name', 'howto']
>>>>>>> remotes/origin/input-web

    def __init__(self, topic_name=None, howto=None):
        self.topic_name = topic_name
        self.howto = howto

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
