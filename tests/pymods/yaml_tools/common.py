'''
Define classes useful in several places and structures required by other modules.
'''
from __future__ import print_function, division, unicode_literals

import re
import sys
import numpy as np

from .abinit_iterators import ITERATOR_RANKS

re_word = re.compile(r'[a-zA-Z0-9_]+')


PY3 = sys.version_info[0] >= 3

if PY3:
    string = str
    basestring = str
else:
    string = unicode
    basestring = basestring


def get_yaml_tag(cls):
    return getattr(cls, '_' + cls.__name__.lstrip('_') + '__yaml_tag', cls.__name__)


def normalize_attr(string):
    return '_'.join(re_word.findall(string))  # .lower()


class BaseDictWrapper(object):
    '''
    Allow attribute access and key access to the values of dictionary to
    keep a consistent behaviour with AutoMap structures. It does not
    inherit from dict but it implements the complete interface.
    '''
    is_dict_like = True

    def __init__(self, d={}, **kwargs):
        for attr in d:
            self[attr] = d[attr]

        for attr in kwargs:
            self[attr] = kwargs[attr]

    def get(self, key, default=None):
        if isinstance(key, basestring):
            key = normalize_attr(key)
        if key in self.__dict__:
            elem = self.__dict__[key]
        else:
            elem = default
        if type(elem) is dict:
            return BaseDictWrapper(elem)
        else:
            return elem

    def __contains__(self, key):
        if isinstance(key, basestring):
            key = normalize_attr(key)
        return key in self.__dict__

    def __getitem__(self, key):
        if isinstance(key, basestring):
            nkey = normalize_attr(key)
        else:
            nkey = key
        if nkey not in self.__dict__:
            raise KeyError(key)
        elem = self.__dict__[nkey]
        if type(elem) is dict:
            return BaseDictWrapper(elem)
        else:
            return elem

    def __setitem__(self, key, val):
        if isinstance(key, basestring):
            key = normalize_attr(key)
        if type(val) is dict:
            val = BaseDictWrapper(val)
        self.__dict__[key] = val

    def __delitem__(self, key):
        nkey = normalize_attr(key)
        if nkey not in self.__dict__:
            raise KeyError(key)
        del self.__dict__[nkey]

    def __repr__(self):
        r = type(self).__name__ + '('
        for attr, val in self.__dict__.items():
            r += '{}={}, '.format(attr, val)
        return r[:-2] + ')'

    def __iter__(self):
        for key in self.__dict__:
            yield key

    def keys(self):
        return self.__dict__.keys()

    def items(self):
        return self.__dict__.items()


class Undef(float):
    '''
        Represent the magic number undef.
    '''
    _is_undef = True
    yaml_pattern = re.compile('undef')

    @staticmethod
    def is_undef(obj):
        return getattr(obj, '_is_undef', False)

    @staticmethod
    def __new__(cls):
        return super(Undef, cls).__new__(cls, 'nan')

    def __eq__(self, other):
        return getattr(other, '_is_undef', False)

    def __repr__(self):
        return 'undef'

    @classmethod
    def from_scalar(cls, scal):
        return cls()

    def to_scalar(self):
        return 'undef'


class FailDetail(object):
    def __init__(self, details):
        self.details = details

    def __bool__(self):
        '''
            As a fail it is always Falsy
        '''
        return False


class BaseArray(np.ndarray):
    '''
        Define a base class for YAML tags converted to numpy compatible
        objects. Can be used for converting any YAML array of number of any
        dimension into a numpy compatible array.
    '''

    # attribute to identify the class without relying on isinstance (unreliable
    # because of sys.path manipulation)
    _is_base_array = True

    # Short tag name
    __yaml_tag = 'Array'

    # by default we want to treat this as a coherent object and do not check
    # values individualy
    has_no_child = True

    def __init__(self, *args, **kwargs):
        # numpy ndarray does not have __init__
        # everything is done in __new__
        self._has_undef = False

    @classmethod
    def from_seq(cls, s):
        def check_undef(s):
            '''
                Look for Undef in the original list because numpy convert it to
                nan
            '''
            if hasattr(s, '__iter__'):
                for el in s:
                    if check_undef(el):
                        return True
                return False
            else:
                return Undef.is_undef(s)

        new = np.array(s).view(cls)
        new._has_undef = check_undef(s)
        return new

    def to_seq(self):
        # conversion have to be explicit because numpy float are not
        # recognised as float by yaml
        def to_list(arr):
            if len(arr.shape) > 1:
                return [to_list(line) for line in arr]
            else:
                return [float(f) for f in arr]
        return to_list(self)


class IterStart(object):
    '''
        Mark the begining of a iteration of a given iterator.
    '''
    # Don't do this at home, trick to workaround the custom sys.path
    _is_iter_start = True

    def __init__(self, iterator, iteration):
        self.iterator = iterator
        self.iteration = iteration

    @classmethod
    def from_map(cls, d):
        iterator = max(d.keys(), key=lambda x: ITERATOR_RANKS[x])
        iteration = d[iterator]
        return cls(iterator, iteration)

    def to_map(self):
        return {self.iterator: self.iteration}

    def __repr__(self):
        return 'IterStart({}={})'.format(self.iterator, self.iteration)
