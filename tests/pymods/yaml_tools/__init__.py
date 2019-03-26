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
    def yaml_parse(content, catch=True, *args, **kwargs):
        from . import structures
        if catch:
            try:
                return yaml.load(content, *args, Loader=yaml.Loader, **kwargs)
            except yaml.YAMLError as e:
                return CorruptedDocument(e, context=content)
        else:
            return yaml.load(content, *args, Loader=yaml.Loader, **kwargs)

    yaml_print = yaml.dump


class CorruptedDocument(object):
    '''
        Replace the YAML parser output when it fail.
    '''
    # trick to workaround the custom sys.path
    _is_corrupted_doc = True

    def __init__(self, error, context=''):
        self.error = error
        self.context = context

    def __repr__(self):
        return "<Corrupted document {}: {}>".format(
            self.error.__class__.__name__, str(self.error))


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
            self._obj = yaml_parse('\n'.join(self.lines))
        else:
            raise NoYAMLSupportError('Try to access YAML document but YAML is'
                                     ' not available in this environment.')
        if isinstance(self._obj, CorruptedDocument):
            self._corrupted = True

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
