'''
Define all error types used in yaml_tools.
All errors must inherit from YAMLTestError.
'''
from __future__ import print_function, division, unicode_literals


class YAMLTestError(Exception):
    """Base class for all other errors."""


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


class AlreadyRegisteredTagError(ConfigParserError):
    def __init__(self, tag):
        msg = 'Attempt to register {} twice.'
        super(ConfigParserError, self).__init__(msg.format(tag))


###############################################################################
class ConfigError(YAMLTestError):
    pass


class ValueTypeError(TypeError, ConfigError):
    def __init__(self, name, exp, found):
        msg = ('The value found in config does not match the type expected for'
               ' {}. Expected {} and found {} of type {}.')
        super(TypeError, self).__init__(msg.format(name, exp, found, type(found)))


class InvalidNodeError(ConfigError):
    def __init__(self, name, value):
        msg = ('The node labeled {} is not a known parameter or constraint and'
               ' have not the form of a specialisation. Value: {}')
        super(InvalidNodeError, self).__init__(msg.format(name, value))


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
        YAMLTestError.__init__(msg)


class UntaggedDocumentError(InputFileError):
    def __init__(self, line):
        msg = ('This document does not have a tag. It cannot be identified.')
        super(UntaggedDocumentError, self).__init__(line, msg)


class TagMismatchError(InputFileError):
    def __init__(self, line, expected, found):
        msg = ('This was supposed to be tagged {} but it was {}.'
               .format(expected, found))
        super(TagMismatchError, self).__init__(line, msg)


class DuplicateDocumentError(InputFileError):
    def __init__(self, line, id):
        msg = ('There are two document with the same tag and iteration'
               ' state ({}). Please change the tag of one of them to make it'
               ' unique.').format(id)
        super(DuplicateDocumentError, self).__init__(line, msg)
