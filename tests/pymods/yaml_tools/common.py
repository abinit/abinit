'''
    Define classes usefull in several places.
'''
from __future__ import print_function, division, unicode_literals
import re
import sys

re_word = re.compile(r'[a-zA-Z0-9_]+')


PY3 = sys.version_info[0] == 3

if PY3:
    string = str
    basestring = str
else:
    string = unicode
    basestring = basestring


def normalize_attr(string):
    return '_'.join(re_word.findall(string))  # .lower()


class BaseDictWrapper(object):
    '''
        Allow attribute access and key access to the values of dictionary to
        keep a consistent behaviour with AutoMap structures. It does not
        inherit from dict but it implement the complete interface.
    '''
    _is_dict_like = True

    def __init__(self, d={}):
        for attr in d:
            self[attr] = d[attr]

    def get(self, key, default=None):
        nkey = normalize_attr(key)
        if nkey in self.__dict__:
            elem = self.__dict__[nkey]
        else:
            elem = default
        if type(elem) is dict:
            return BaseDictWrapper(elem)
        else:
            return elem

    def __contains__(self, key):
        nkey = normalize_attr(key)
        return nkey in self.__dict__

    def __getitem__(self, key):
        nkey = normalize_attr(key)
        if nkey not in self.__dict__:
            raise KeyError(key)
        elem = self.__dict__[nkey]
        if type(elem) is dict:
            return BaseDictWrapper(elem)
        else:
            return elem

    def __setitem__(self, key, val):
        if type(val) is dict:
            val = BaseDictWrapper(val)
        self.__dict__[normalize_attr(key)] = val

    def __delitem__(self, key):
        nkey = normalize_attr(key)
        if nkey not in self.__dict__:
            raise KeyError(key)
        del self.__dict__[nkey]

    def __repr__(self):
        r = self.__class__.__name__ + '('
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
