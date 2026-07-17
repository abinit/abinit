"""
Define all error types used in yaml_tools.
All errors must inherit from YAMLTestError.
"""


class YAMLTestError(Exception):
    """Base class for all other errors."""


class ConfigContextError(YAMLTestError):
    def __init__(self, path):
        """
        Initialize ConfigContextError.

        Args:
            path (list): The path in the config tree where the error occurred.
        """
        spath = ".".join(path)
        msg = (f"Tried to enter a None context in the config tree at {spath}"
               )
        super().__init__(msg)


class NoYAMLSupportError(YAMLTestError):
    """Raised when Yaml library is not installed."""


###############################################################################
class ConfigParserError(YAMLTestError):
    pass


class UnknownParamError(ConfigParserError):
    def __init__(self, cons, param):
        """
        Initialize UnknownParamError.

        Args:
            cons (str): Name of the constraint.
            param (str): Name of the unknown parameter.
        """
        msg = ('Encounterd an unknown parameter name "{}"'
               ' when registering constraint "{}".')
        super().__init__(msg.format(param, cons))


class AlreadyRegisteredTagError(ConfigParserError):
    def __init__(self, tag):
        """
        Initialize AlreadyRegisteredTagError.

        Args:
            tag (str): The tag that was already registered.
        """
        msg = "Attempt to register {} twice."
        super(ConfigParserError, self).__init__(msg.format(tag))


###############################################################################
class ConfigError(YAMLTestError):
    pass


class ValueTypeError(TypeError, ConfigError):
    def __init__(self, name, exp, found):
        """
        Initialize ValueTypeError.

        Args:
            name (str): Name of the parameter or constraint.
            exp (type): The expected type.
            found: The value found in the configuration.
        """
        msg = ("The value found in config does not match the type expected for"
               " {}. Expected {} and found {} of type {}.")
        super(TypeError, self).__init__(msg.format(name, exp, found, type(found)))


class InvalidNodeError(ConfigError):
    def __init__(self, name, value):
        """
        Initialize InvalidNodeError.

        Args:
            name (str): The label of the invalid node.
            value: The value of the invalid node.
        """
        msg = ("The node labeled {} is not a known parameter or constraint and"
               " have not the form of a specialisation. Value: {}")
        super().__init__(msg.format(name, value))


class EmptySetError(ConfigError):
    def __init__(self, obj):
        """
        Initialize EmptySetError.

        Args:
            obj: The object used to attempt creating an empty set.
        """
        msg = "User tried to create an empty set with {}."
        super().__init__(msg.format(obj))


class NotOrderedOverlappingSetError(ConfigError):
    def __init__(self, set1, set2):
        msg = "{} and {} are overlapping but cannot be ordered."
        super().__init__(msg.format(set1,
                                                                       set2))


class IllegalFilterNameError(ConfigError):
    def __init__(self, name):
        msg = "{} is a reserved name, you cannot use it as a filter name."
        super().__init__(msg.format(name))


class MissingCallbackError(ConfigError):
    def __init__(self, obj, method):
        msg = f"{obj} does not expose a {method} method."
        super().__init__(msg)


###############################################################################
class InputFileError(YAMLTestError):
    def __init__(self, line, msg):
        msg = f"In input file at line {line}:\n{msg}"
        super().__init__(self, msg)


class NoIteratorDefinedError(InputFileError):
    def __init__(self, doc):
        msg = (f"No iterator have been found before the first document {doc.obj}."
               )
        super().__init__(doc.start + 1, msg)


class NotAvailableTagError(InputFileError):
    def __init__(self, tag, msg):
        msg = (f"Tag {tag} is not available in this installation : {msg}"
               )
        YAMLTestError.__init__(msg)


class UntaggedDocumentError(InputFileError):
    def __init__(self, line):
        msg = ("This document does not have a tag. It cannot be identified.")
        super().__init__(line, msg)


class TagMismatchError(InputFileError):
    def __init__(self, line, expected, found):
        msg = (f"This was supposed to be tagged {expected} but it was {found}."
               )
        super().__init__(line, msg)


class DuplicateDocumentError(InputFileError):
    def __init__(self, line, id):
        msg = ("There are two document with the same tag and iteration"
               f" state ({id}). Please change the tag of one of them to make it"
               " unique.")
        super().__init__(line, msg)
