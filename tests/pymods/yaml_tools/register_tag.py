"""
This module provides facilities to use YAML formatted data.
It defines several decorators to easily create YAML compatible
classes which are used both when parsing YAML formatted data
and when writing YAML formatted data.
"""

import re
import warnings
from inspect import ismethod

import yaml

from . import Loader
from .common import BaseDictWrapper, get_yaml_tag
from .errors import AlreadyRegisteredTagError, NotAvailableTagError

known_tags = set()


def reserve_tag(tag):
    """
    Prevent multiple registration of the same tag.

    Args:
        tag (str): The tag to reserve.

    Raises:
        AlreadyRegisteredTagError: If the tag is already in `known_tags`.
    """
    if tag in known_tags:
        raise AlreadyRegisteredTagError(tag)
    known_tags.add(tag)


def yaml_map(cls):
    """
    Register a class as a YAML mapping (!Tag).

    The class must expose `from_map` and `to_map`.
    `from_map` returns a valid instance and can be a classmethod or instance method.

    Args:
        cls (type): The class to register.

    Returns:
        type: The registered class.
    """
    tag = "!" + get_yaml_tag(cls)

    reserve_tag(tag)

    def constructor(loader, node):
        map = dict(loader.construct_mapping(node, deep=True))
        if ismethod(cls.from_map):
            return cls.from_map(map)
        return cls().from_map(map)

    def representer(dumper, data):
        return dumper.represent_mapping(tag, data.to_map())

    yaml.add_constructor(tag, constructor, Loader=Loader)
    yaml.add_representer(cls, representer)

    return cls


def yaml_seq(cls):
    """
    Register a class as a YAML sequence (!Tag).

    The class must expose `from_seq` and `to_seq`.
    `from_seq` returns a valid instance and can be a classmethod or instance method.

    Args:
        cls (type): The class to register.

    Returns:
        type: The registered class.
    """
    tag = "!" + get_yaml_tag(cls)

    reserve_tag(tag)

    def constructor(loader, node):
        seq = list(loader.construct_sequence(node, deep=True))
        if ismethod(cls.from_seq):
            return cls.from_seq(seq)
        return cls().from_seq(seq)

    def representer(dumper, data):
        return dumper.represent_sequence(tag, data.to_seq())

    yaml.add_constructor(tag, constructor, Loader=Loader)
    yaml.add_representer(cls, representer)

    return cls


def yaml_scalar(cls):
    """
    Register a class as a YAML scalar (!Tag).

    The class must expose `from_scalar` and `to_scalar`.
    `from_scalar` returns a valid instance and can be a classmethod or instance method.

    Args:
        cls (type): The class to register.

    Returns:
        type: The registered class.
    """
    tag = "!" + get_yaml_tag(cls)

    reserve_tag(tag)

    def constructor(loader, node):
        scalar = loader.construct_scalar(node)
        if ismethod(cls.from_scalar):
            return cls.from_scalar(scalar)
        return cls().from_scalar(scalar)

    def representer(dumper, data):
        return dumper.represent_scalar(tag, data.to_scalar())

    yaml.add_constructor(tag, constructor, Loader=Loader)
    yaml.add_representer(cls, representer)

    return cls


def auto_map(Cls):
    """
    Automatically append `from_map`, `to_map`, and `__repr__` to a class.

    Requires `__getitem__` and `__setitem__` to be defined. Attribute names
    are normalized for property-like access.

    Args:
        Cls (type): The class to enhance.

    Returns:
        type: The enhanced AutoMap class.

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
    ... a['attr .w. --spaces']
    82
    """
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
    """
    Register a class as a YAML mapping with auto-generated map methods.

    Args:
        cls (type): The class to register.

    Returns:
        type: The registered class.
    """
    return yaml_map(auto_map(cls))


def yaml_implicit_scalar(cls):
    """
    Register a class as a YAML scalar with an implicit pattern.

    The class must expose `from_scalar` and `to_scalar`, and have a
    `yaml_pattern` attribute (string or regex).

    Args:
        cls (type): The class to register.

    Returns:
        type: The registered class.
    """
    yaml_scalar(cls)  # register the constructor and the representer
    tag = "!" + get_yaml_tag(cls)

    re_pattern = cls.yaml_pattern
    if not hasattr(re_pattern, "match"):
        re_pattern = re.compile(re_pattern)
    # register the implicit pattern
    yaml.add_implicit_resolver(tag, re_pattern, Loader=Loader)
    return cls


def yaml_not_available_tag(tag, reason, fatal=False):
    """
    Register tag with a given tag but trigger a warning if fatal == false
    or an error if fatal == True. Use `reason` as the message
    to give more detail to the user. If fatal is False then the returned object
    is an empty dictionary.
    """
    msg = f"The tag !{tag} is used but is not available:\n{reason}"

    reserve_tag("!" + tag)

    def constructor(loader, node):
        if fatal:
            raise NotAvailableTagError(msg)
        warnings.warn(msg)
        return {"_not_available": True}
    yaml.add_constructor("!" + tag, constructor, Loader=Loader)
