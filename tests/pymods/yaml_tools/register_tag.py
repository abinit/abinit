'''
    This module provide some facilities to use YAML formated data.
    It defines several decorators to easily create YAML compatible
    classes witch will be used both in parsing YAML formated data
    and in writing YAML formated data.
'''
import re
from inspect import ismethod

try:
    import yaml
    is_available = True
except ImportError:
    is_available = False

if is_available:
    def yaml_map(Cls):
        '''
            Register a class as a known tag for YAML.
            The class must expose methods from_map and to_map.
            from_map must return a valid instance of Cls.
            It can be a classmethod, a staticmethod or a normal
            method but in the last case 'Cls()' must be a valid
            initialisation.
        '''
        tag = '!' + getattr(Cls, Cls.__name__ + '_yaml_tag', Cls.__name__)

        def constructor(loader, node):
            map = dict(loader.construct_mapping(node, deep=True))
            if ismethod(Cls.from_map):
                return Cls.from_map(map)
            else:
                return Cls().from_map(map)

        def representer(dumper, data):
            return dumper.represent_mapping(tag, data.to_map())

        yaml.add_constructor(tag, constructor)
        yaml.add_representer(Cls, representer)

        return Cls

    def yaml_seq(Cls):
        '''
            Register a class as a known tag for YAML.
            The class must expose methods from_seq and to_seq.
            from_seq must return a valid instance of Cls.
            It can be a classmethod, a staticmethod or a normal
            method but in the last case 'Cls()' must be a valid
            initialisation.
        '''
        tag = '!' + getattr(Cls, Cls.__name__ + '_yaml_tag', Cls.__name__)

        def constructor(loader, node):
            seq = list(loader.construct_sequence(node, deep=True))
            if ismethod(Cls.from_seq):
                return Cls.from_seq(seq)
            else:
                return Cls().from_seq(seq)

        def representer(dumper, data):
            return dumper.represent_sequence(tag, data.to_seq())

        yaml.add_constructor(tag, constructor)
        yaml.add_representer(Cls, representer)

        return Cls

    def yaml_scalar(Cls):
        '''
            Register a class as a known tag for YAML.
            The class must expose methods from_scalar to_scalar.
            from_scalar must return a valid instance of Cls.
            It can be a classmethod, a staticmethod or a normal
            method but in the last case 'Cls()' must be a valid
            initialisation.
        '''
        tag = '!' + getattr(Cls, Cls.__name__ + '_yaml_tag', Cls.__name__)

        def constructor(loader, node):
            scalar = loader.construct_scalar(node)
            if ismethod(Cls.from_scalar):
                return Cls.from_scalar(scalar)
            else:
                return Cls().from_seq(scalar)

        def representer(dumper, data):
            return dumper.represent_scalar(tag, data.to_scalar())

        yaml.add_constructor(tag, constructor)
        yaml.add_representer(Cls, representer)

        return Cls

    def auto_map(Cls):
        '''
            Automatically append methods from_map, to_map and __repr__ to a
            class intended to be used with YAML tag
        '''
        class AutoMap(Cls):
            @classmethod
            def from_map(cls, d):
                new = cls()
                for attr in new.__dict__:
                    if attr in d:
                        new.__dict__[attr] = d[attr]
                return new

            def to_map(self):
                return self.__dict__

            def __repr__(self):
                r = Cls.__name__ + '('
                for attr, val in self.__dict__.items():
                    r += '{}={}, '.format(attr, val)
                return r[:-2] + ')'

        AutoMap.__name__ = Cls.__name__
        return AutoMap

    def yaml_auto_map(Cls):
        '''
                @yaml_auto_map
            is equivalent to:
                @yaml_map
                @auto_map
        '''
        return yaml_map(auto_map(Cls))

    def yaml_implicit_scalar(Cls):
        '''
            Register a class as a known tag for YAML and a pattern
            wich imply this tag.
            The class must expose methods from_scalar to_scalar.
            Moreover it must have a class attribute yaml_pattern matching
            all recognised strings representing this kind of data.
        '''
        yaml_scalar(Cls)  # register the constructor and the representer
        tag = '!' + getattr(Cls, Cls.__name__ + '_yaml_tag', Cls.__name__)

        re_pattern = Cls.yaml_pattern
        if not hasattr(re_pattern, 'match'):
            re_pattern = re.compile(re_pattern)
        yaml.add_implicit_resolver(tag, re_pattern)  # register the mplicit pattern
        return Cls
