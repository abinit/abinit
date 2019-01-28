from yaml import load, YAMLObject


def yamlobject(ytag='', **attributes):
    '''
        Generate a suitable object for using tags in yaml source.
        Providing named arguments associated with default values
        your will create a default constructor and a default
        __repr__ method for your YAML tag. By default the YAML tag
        is the name of the class you defined.
        If needed you replace this tag with the optional argument
        ytag (do not prepend the '!').
        If your class __repr__ method is not inherited from object
        it will not be overwritten.
        Do not bypass this mecanism for the definition of yaml_tag
        or the pyyaml magic will mess everything.
    '''
    def dec(cls):
        assert issubclass(cls, object), 'yaml_object class sould at least inherit from object'

        class YamlObj(YAMLObject, cls):
            if ytag == '':
                yaml_tag = '!' + cls.__name__
            else:
                yaml_tag = '!' + ytag

            @classmethod
            def __new__(clsb, clsc):
                '''
                Setup the default values directly in __new__ because yaml bypass __init__
                '''
                self = cls.__new__(clsb)
                for attr, default in attributes.items():
                    self.__setattr__(attr, default)
                return self

        if cls.__repr__ == object.__repr__:
            def repr(self):
                r = cls.__name__ + '('
                for attr in attributes:
                    r += '{}={}, '.format(attr, self.__getattribute__(attr))
                r = r[:-2] + ')'
                return r
            YamlObj.__repr__ = repr
        return YamlObj
    return dec


yaml_parse = load
