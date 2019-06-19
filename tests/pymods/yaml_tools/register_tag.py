'''
    This module provides facilities to use YAML formatted data.
    It defines several decorators to easily create YAML compatible
    classes which will be used both when parsing YAML formatted data
    and whwn writing YAML formatted data.
'''
from __future__ import print_function, division, unicode_literals

import re
import warnings
import yaml

from inspect import ismethod
from . import Loader
from .common import BaseDictWrapper
from .errors import NotAvailableTagError


def yaml_tag_mangle(cls):
    """Return the mangled name of the attribute __yaml_tag."""
    return '_' + cls.__name__.lstrip('_') + '__yaml_tag'


def yaml_map(cls):
    '''
    Register a class with a given tag in the YAML library.
    The class must expose the methods `from_map` and `to_map`.
    `from_map` must return a valid instance of cls.
    It can be a classmethod or a normal method but in
    the latter case 'cls()' must be a valid initialisation.
    '''
    tag = '!' + getattr(cls, yaml_tag_mangle(cls), cls.__name__)

    def constructor(loader, node):
        map = dict(loader.construct_mapping(node, deep=True))
        if ismethod(cls.from_map):
            return cls.from_map(map)
        else:
            return cls().from_map(map)

    def representer(dumper, data):
        return dumper.represent_mapping(tag, data.to_map())

    yaml.add_constructor(tag, constructor, Loader=Loader)
    yaml.add_representer(cls, representer)

    return cls


def yaml_seq(cls):
    '''
    Register a class with a given tag in the YAML library.
    The class must expose the methods `from_seq` and `to_seq`.
    `from_seq` must return a valid instance of cls.
    It can be a class method or a normal method but in
    the latter case 'cls()' must be a valid initialisation.
    '''
    tag = '!' + getattr(cls, yaml_tag_mangle(cls), cls.__name__)

    def constructor(loader, node):
        seq = list(loader.construct_sequence(node, deep=True))
        if ismethod(cls.from_seq):
            return cls.from_seq(seq)
        else:
            return cls().from_seq(seq)

    def representer(dumper, data):
        return dumper.represent_sequence(tag, data.to_seq())

    yaml.add_constructor(tag, constructor, Loader=Loader)
    yaml.add_representer(cls, representer)

    return cls


def yaml_scalar(cls):
    '''
    Register a class with a given tag in the YAML library.
    The class must expose the methods `from_scalar` and `to_scalar`.
    `from_scalar` must return a valid instance of cls.
    It can be a class method or a normal method but in
    the latter case 'cls()' must be a valid initialisation.
    '''
    tag = '!' + getattr(cls, yaml_tag_mangle(cls), cls.__name__)

    def constructor(loader, node):
        scalar = loader.construct_scalar(node)
        if ismethod(cls.from_scalar):
            return cls.from_scalar(scalar)
        else:
            return cls().from_scalar(scalar)

    def representer(dumper, data):
        return dumper.represent_scalar(tag, data.to_scalar())

    yaml.add_constructor(tag, constructor, Loader=Loader)
    yaml.add_representer(cls, representer)

    return cls


def auto_map(Cls):
    '''
        Automatically append methods from_map, to_map and __repr__ to a
        class intended to be used with YAML tag, provided __getitem__
        and __setitem__ are defined. Attribute names are normalized to 
        be accessible as regular property even if the
        original name contained spaces or special characters. The original
        name can still be used in dict like access.
        Example:

        >>> @auto_map
        ... class A(object):
        ...     pass
        ...
        >>> a = A()
        >>> a['attr w/ spaces'] = 78
        >>> a.attr_w_spaces
        78
        >>> a.attr_w_spaces = 82
        >>> a['attr w/ spaces']
        82
        >>> # be careful, simple normalization imply collisions
        >>> a['attr .w. --spaces']
        82
    '''
    class AutoMap(Cls, BaseDictWrapper):
        @classmethod
        def from_map(cls, d):
            new = cls()
            for attr in d:
                new[attr] = d[attr]
            return new

        def to_map(self):
            return self.__dict__

    AutoMap.__name__ = Cls.__name__
    return AutoMap


def yaml_auto_map(cls):
    '''
    @yaml_auto_map

    is equivalent to:

    @yaml_map
    @auto_map
    '''
    return yaml_map(auto_map(cls))


def yaml_implicit_scalar(cls):
    '''
    Register a class with a given tag in the YAML library and a pattern
    which imply this tag.
    The class must expose methods `from_scalar` and `to_scalar`.
    Moreover it must have a class attribute `yaml_pattern` that
    can be either a string with the regex pattern matching
    all string representing this kind of data or the compiled
    object of this same regex.
    '''
    yaml_scalar(cls)  # register the constructor and the representer
    tag = '!' + getattr(cls, yaml_tag_mangle(cls), cls.__name__)

    re_pattern = cls.yaml_pattern
    if not hasattr(re_pattern, 'match'):
        re_pattern = re.compile(re_pattern)
    # register the implicit pattern
    yaml.add_implicit_resolver(tag, re_pattern, Loader=Loader)
    return cls


def yaml_not_available_tag(tag, reason, fatal=False):
    '''
    Register tag with a given tag but trigger a warning if fatal == false
    or an error if fatal == True. Use `reason` as the message
    to give more detail to the user. If fatal is False then the returned object
    is an empty dictionary.
    '''
    msg = 'The tag !{} is used but is not available:\n{}'.format(tag, reason)

    def constructor(loader, node):
        if fatal:
            raise NotAvailableTagError(msg)
        else:
            warnings.warn(msg)
            return {'_not_available': True}
    yaml.add_constructor('!' + tag, constructor, Loader=Loader)
