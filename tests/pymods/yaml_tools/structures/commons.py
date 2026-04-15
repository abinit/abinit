"""
Define basic structures.
"""
from ..register_tag import yaml_auto_map, yaml_implicit_scalar


@yaml_auto_map
class GenericMap:
    """A generic tag definition for test and example."""


@yaml_implicit_scalar
class YAMLComplex(complex):
    #                 >             [1]                       <
    yaml_pattern = (r"[+-]?(\d+(\.\d*)?|\.\d+)([eEdD][+-]?\d+)?"
                    r" *[+-] *[+-]?(\d+(\.\d*)?|\.\d+)([eEdD][+-]?\d+)?i")
    #                 >  [2] <>                    [3]                <
    # [1] and [3] float with optional sign and exponential notation, will
    # also match integers and .1 like (fortran does not produce this though)
    # [2] + or - with optional blanks around

    @staticmethod
    def __new__(*args, **kwargs):
        return complex.__new__(*args, **kwargs)

    @classmethod
    def from_scalar(cls, scal):
        """
        Convert a YAML complex number scalar to a Python complex object.

        Args:
            scal (str): The scalar string representing the complex number.

        Returns:
            YAMLComplex: The parsed complex number.
        """
        return cls(scal
                   # python always uses double and only recognise E and e
                   .replace("d", "e")
                   .replace("D", "e")
                   # python uses j instead of i (as in electro magnetism)
                   .replace("i", "j")
                   # spaces have to be stripped out around the central + or -
                   .replace(" ", "")
                   # python expects only one + or - in string form
                   .replace("+-", "-")
                   .replace("-+", "-"))

    def to_scalar(self):
        return repr(self)[1:-1]  # remove parentheses


class AbinitMessage:
    _is_abinit_message = True


@yaml_auto_map
class AbinitError(AbinitMessage):
    """Base class for Abinit messages."""
    __yaml_tag = "ERROR"


@yaml_auto_map
class AbinitWarning(AbinitMessage):
    __yaml_tag = "WARNING"


# MG: Is this uses somewhere?
@yaml_auto_map
class AbinitInfo(AbinitMessage):
    __yaml_tag = "INFO"


@yaml_auto_map
class AbinitComment(AbinitMessage):
    __yaml_tag = "COMMENT"


@yaml_auto_map
class DatasetInfo:
    __yaml_tag = "DatasetInfo"


@yaml_auto_map
class BeginCycle:
    __yaml_tag = "BeginCycle"
