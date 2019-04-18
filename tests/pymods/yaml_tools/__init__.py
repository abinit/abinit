'''
    This package gather all tools used by the Abinit test suite for
    manipulating YAML formated data.
'''
from __future__ import print_function, division, unicode_literals

import warnings
from .errors import NoYAMLSupportError

try:
    import yaml
    import numpy  # numpy is also required

    is_available = True

except ImportError:
    warnings.warn('Cannot import numpy or yaml package. Use `pip install numpy'
                  ' pyyaml --user` to install the packages in user mode.')
    is_available = False

try:
    import pandas
    has_pandas = True
except ImportError:
    has_pandas = False
    warnings.warn('Cannot import pandas package. Use `pip install pandas'
                  ' --user` to install the package in user mode.')


if is_available:
    if hasattr(yaml, 'CSafeLoader'):  # use the C binding (faster) if possible
        Loader = yaml.CSafeLoader
    else:
        Loader = yaml.SafeLoader

    from .common import Undef, IterStart

    def yaml_parse(content, *args, **kwargs):
        from .register_tag import yaml_implicit_scalar, yaml_map
        yaml_implicit_scalar(Undef)
        yaml_map(IterStart)

        from . import structures
        return yaml.load(content, *args, Loader=Loader, **kwargs)

    yaml_print = yaml.dump


class Document(object):
    '''
        Represent a document with all its metadata from the original file.
    '''

    def __init__(self, iterators, start, lines):
        self.iterators = iterators
        self.start = start
        self.end = -1
        self.lines = lines
        self._obj = None
        self._corrupted = False

    def _parse(self):
        if is_available:
            content = '\n'.join(self.lines)
            try:
                self._obj = yaml_parse(content)
            except yaml.YAMLError as e:
                self._obj = e
                self._corrupted = True
        else:
            raise NoYAMLSupportError('Try to access YAML document but YAML is'
                                     ' not available in this environment.')

    @property
    def obj(self):
        if self._obj is None:
            self._parse()
        return self._obj

    @property
    def corrupted(self):
        if self._obj is None:
            self._parse()
        return self._corrupted
