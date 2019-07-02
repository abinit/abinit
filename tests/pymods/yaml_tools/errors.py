from __future__ import print_function, division, unicode_literals


class YAMLTestError(Exception):
    """Base class."""


class ConfigContextError(YAMLTestError):
    def __init__(self, path):
        spath = '.'.join(path)
        msg = ('Tried to enter a None context in the config tree at {}'
               .format(spath))
        super(ConfigContextError, self).__init__(msg)


class NoYAMLSupportError(YAMLTestError):
    """Raised when Yaml library is not installed."""


###############################################################################
class ConfigParserError(YAMLTestError):
    pass


class UnknownParamError(ConfigParserError):
    def __init__(self, cons, param):
        msg = ('Encounterd an unknown parameter name "{}"'
               ' when registering constraint "{}".')
        super(UnknownParamError, self).__init__(msg.format(param, cons))


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
        super(InputFileError, self).__init__(msg.format(name, value))


class EmptySetError(ConfigError):
    def __init__(self, obj):
        msg = 'User tried to create an empty set with {}.'
        super(EmptySetError, self).__init__(msg.format(obj))


class NotOrderedOverlappingSetError(ConfigError):
    def __init__(self, set1, set2):
        msg = '{} and {} are overlapping but cannot be ordered.'
        super(NotOrderedOverlappingSetError, self).__init__(msg.format(set1,
                                                                       set2))


class IllegalFilterNameError(ConfigError):
    def __init__(self, name):
        msg = '{} is a reserved name, you cannot use it as a filter name.'
        super(IllegalFilterNameError, self).__init__(msg.format(name))


class MissingCallbackError(ConfigError):
    def __init__(self, obj, method):
        msg = '{} does not expose a {} method.'.format(obj, method)
        super(MissingCallbackError, self).__init__(msg)


###############################################################################
class InputFileError(YAMLTestError):
    def __init__(self, line, msg):
        msg = 'In input file at line {}:\n{}'.format(line, msg)
        super(InputFileError, self).__init__(self, msg)


class NoIteratorDefinedError(InputFileError):
    def __init__(self, doc):
        msg = ('No iterator have been found before the first document {}.'
               .format(doc.obj))
        super(NoIteratorDefinedError, self).__init__(doc.start + 1, msg)


class NotAvailableTagError(InputFileError):
    def __init__(self, tag, msg):
        msg = ('Tag {} is not available in this installation : {}'
               .format(tag, msg))
        YAMLTestError.__init__(self, msg)


class UnlabeledDocumentError(InputFileError):
    def __init__(self, line):
        msg = ('This document does not have a label field. It cannot be'
               ' identified.')
        InputFileError.__init__(self, line, msg)


class DuplicateDocumentError(InputFileError):
    def __init__(self, line, id):
        msg = ('There are two document with the same label and iteration'
               ' state ({}). Please change the label of one of them to make it'
               ' unique.').format(id)
        InputFileError.__init__(self, line, msg)
