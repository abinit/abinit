'''
    Define basic structures.
'''
from __future__ import print_function, division, unicode_literals
from ..register_tag import yaml_auto_map, yaml_implicit_scalar


@yaml_implicit_scalar
class YAMLComplex(complex):
    #                 >             [1]                       <
    yaml_pattern = (r'[+-]?(\d+(\.\d*)?|\.\d+)([eEdD][+-]?\d+)?'
                    r' *[+-] *[+-]?(\d+(\.\d*)?|\.\d+)([eEdD][+-]?\d+)?i')
    #                 >  [2] <>                    [3]                <
    # [1] and [3] float with optional sign and exponential notation, will
    # also match integers and .1 like (fortran does not produce this though)
    # [2] + or - with optional blanks around

    @staticmethod
    def __new__(*args, **kwargs):
        return complex.__new__(*args, **kwargs)

    @classmethod
    def from_scalar(cls, scal):
        return cls(scal
                   # python always use double and only recognise E and e
                   .replace('d', 'e')
                   .replace('D', 'e')
                   # python use j instead of i (as in electro magnetism)
                   .replace('i', 'j')
                   # spaces have to be striped around the central + or -
                   .replace(' ', '')
                   # python expect only on + or - in string form
                   .replace('+-', '-')
                   .replace('-+', '-'))

    def to_scalar(self):
        return repr(self)[1:-1]  # remove paranthesis


class AbinitMessage(object):
    _is_abinit_message = True


@yaml_auto_map
class AbinitError(AbinitMessage):
    __yaml_tag = 'ERROR'


@yaml_auto_map
class AbinitWarning(AbinitMessage):
    __yaml_tag = 'WARNING'


@yaml_auto_map
class AbinitInfo(object):
    __yaml_tag = 'INFO'


@yaml_auto_map
class AbinitComment(AbinitMessage):
    __yaml_tag = 'COMMENT'
