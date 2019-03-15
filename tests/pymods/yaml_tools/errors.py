from __future__ import print_function, division, unicode_literals


class YAMLTestError(Exception):
    pass


class ConfigContextError(YAMLTestError):
    def __init__(self, path):
        spath = '.'.join(path)
        msg = ('Tried to enter a None context in the config tree at {}'
               .format(spath))
        super(ConfigContextError, self).__init__(self, msg)


###############################################################################
class ConfigParserError(YAMLTestError):
    pass


class UnknownParamError(ConfigParserError):
    def __init__(self, cons, param):
        msg = ('Encounterd an unknown parameter name "{}"'
               ' when registering constraint "{}".')
        super(UnknownParamError, self).__init__(self, msg.format(param, cons))


###############################################################################
class ConfigError(YAMLTestError):
    pass


class ValueTypeError(TypeError, ConfigError):
    def __init__(self, name, exp, found):
        msg = ('The value found in config does not match the type expected for'
               ' {}. Expected {} and found {} of type {}.')
        TypeError.__init__(self, msg.format(name, exp, found, type(found)))


class InvalidNodeError(ConfigError):
    def __init__(self, name, value):
        msg = ('The node labeled {} is not a known parameter or constraint and'
               ' have not the form of a specialisation. Value: {}')
        ConfigError.__init__(self, msg.format(name, value))


class AlreadySetKeyError(ConfigError):
    def __init__(self, key):
        msg = 'The {} key have already been set in this tree.'
        ConfigError.__init__(self, msg.format(key))


class EmptySetError(ConfigError):
    def __init__(self, obj):
        msg = 'User tried to create an empty set with {}.'
        ConfigError.__init__(self, msg.format(obj))


class NotOrderedOverlappingSetError(ConfigError):
    def __init__(self, set1, set2):
        msg = '{} and {} are overlapping but cannot be ordered.'
        ConfigError.__init__(self, msg.format(set1, set2))


class IllegalFilterNameError(ConfigError):
    def __init__(self, name):
        msg = '{} is a reserved name, you cannot use it as a filter name.'
        ConfigError.__init__(self, msg.format(name))


###############################################################################
class InputFileError(YAMLTestError):
    def __init__(self, line, msg):
        msg = 'In input file at line {}:\n{}'.format(line, msg)
        super(InputFileError, self).__init__(self, msg)


class NoIteratorDefinedError(InputFileError):
    def __init__(self, doc):
        msg = ('No iterator have been found before the first document {}.'
               .format(doc['obj']))
        super(InputFileError, self).__init__(self, doc['start']+1, msg)
