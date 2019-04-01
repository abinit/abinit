'''
    Define some structures to be found in Abinit YAML formated output extending
    the possible operations on the extracted data. All structures that have
    children have to hinerit from BaseDataStructure
'''
from __future__ import print_function, division, unicode_literals
import numpy as np
from . import has_pandas
from .common import FailDetail
from .abinit_iterators import ITERATOR_RANKS
from .register_tag import (
    yaml_map, yaml_seq, yaml_auto_map, yaml_implicit_scalar, yaml_scalar,
    yaml_not_available_tag
)


@yaml_map
class IterStart(object):
    '''
        Mark the begining of a iteration of a given iterator.
    '''
    # Don't do this at home, trick to workaround the custom sys.path
    _is_iter_start = True

    def __init__(self, iterator, iteration):
        self.iterator = iterator
        self.iteration = iteration

    @classmethod
    def from_map(cls, d):
        iterator = max(d.keys(), key=lambda x: ITERATOR_RANKS[x])
        iteration = d[iterator]
        return cls(iterator, iteration)

    def to_map(self):
        return {self.iterator: self.iteration}

    def __repr__(self):
        return 'IterStart({}={})'.format(self.iterator, self.iteration)


@yaml_auto_map
class Etot(object):
    __yaml_tag = 'ETOT'

    def __init__(self, label='nothing', comment='no comment'):
        self.label = label
        self.comment = comment

    @classmethod
    def from_map(cls, map):
        new = super(Etot, cls).from_map(map)
        new.components = {
            name: value for name, value in new.__dict__.items()
            if name not in [
                'Etotal',
                'label',
                'comment',
                'Band energy',
                'Total energy(eV)'
            ]
        }
        return new


@yaml_seq
class AutoNumpy(np.ndarray):
    '''
        Define a base class for YAML tags converted to numpy compatible
        objects.  Can be used for converting any YAML array of number of any
        dimension into a numpy compatible array.
    '''
    __yaml_tag = 'Array'

    # by default we want to treat this as a coherent object and do not check
    # values individualy
    _has_no_child = True

    @classmethod
    def from_seq(cls, s):
        return np.array(s).view(cls)

    def to_seq(self):
        # conversion have to be explicit because numpy float are not
        # recognised as float by yaml
        def to_list(arr):
            if len(arr.shape) > 1:
                return [to_list(line) for line in arr]
            else:
                return [float(f) for f in arr]
        return to_list(self)


@yaml_seq
class Atoms3D(AutoNumpy):
    pass


@yaml_seq
class CartForces(Atoms3D):
    pass


@yaml_seq
class Matrix33(AutoNumpy):
    '''
        Define a matrix of shape (3, 3) compatible with numpy arrays and
        with YAML tags.
    '''
    def __init__(self, shape=(3, 3), *args, **kwargs):
        # numpy ndarray does not have __init__
        # everything is done in __new__
        assert shape == (3, 3)

    @classmethod
    def from_seq(cls, s):
        new = super(Matrix33, cls).from_seq(s)
        assert new.shape == (3, 3)
        return new


@yaml_seq
class Tensor(Matrix33):
    def is_symetric(self):
        for i in range(3):
            for j in range(i, 3):
                if self[i, j] != self[j, i]:
                    return False
        return True

    def is_anti_symetric(self):
        for i in range(3):
            for j in range(i, 3):
                if self[i, j] != -self[j, i]:
                    return False
        return True


@yaml_seq
class Tensor32(Matrix33):
    '''
        Define a matrix of shape (3, 3) compatible with numpy arrays
        and using a Voight notation in YAML.
    '''
    @classmethod
    def from_seq(cls, s):
        assert len(s) == 2
        assert len(s[0]) == 3
        assert len(s[1]) == 3
        arr = np.zeros((3, 3))
        arr[0, 0] = s[0][0]
        arr[1, 1] = s[0][1]
        arr[2, 2] = s[0][2]

        arr[1, 2] = arr[2, 1] = s[1][0]
        arr[0, 2] = arr[2, 0] = s[1][1]
        arr[0, 1] = arr[1, 0] = s[1][2]
        return arr.view(cls)

    def to_seq(self):
        seq = [[0.0] * 3, [0.0] * 3]
        seq[0][0] = float(self[0, 0])  # we have to explicitly convert values
        seq[0][1] = float(self[1, 1])  # to float for yaml to recognise it
        seq[0][2] = float(self[2, 2])  # properly because those values are

        seq[1][0] = float(self[1, 2])  # float compatible but of a numpy
        seq[1][1] = float(self[0, 2])  # custom type
        seq[1][2] = float(self[0, 1])
        return seq


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


@yaml_auto_map
class GwSigma(object):
    pass


if has_pandas:
    from pandas import read_csv, DataFrame
    from pandas.compat import StringIO

    @yaml_scalar
    class Table(DataFrame):
        _is_dict_like = True
        table_sep = r'\s+'

        @classmethod
        def from_scalar(cls, scal):
            return cls(read_csv(StringIO(scal), sep=cls.table_sep))

        def to_scalar(self):
            return self.to_string()

    @yaml_scalar
    class GwSigmaData(Table):
        pass

    @yaml_scalar
    class EtotIters(Table):
        _is_dict_like = False  # prevent tester from going inside by itself

        residues = {
            'deltaE(h)',
            'residm',
            'vres2'
        }

        def last_iter(self, other, **opts):
            tol = opts.get('tol', 1.0e-10)
            ceil = opts.get('ceil', 1.0e-10)
            tol_iter = opts.get('tol_iter', 5)

            def chk_tol(a, b):
                return abs(a - b) / (abs(a) + abs(b)) < tol

            def chk_ceil(b):
                return abs(b) < ceil

            for key in self:
                oserie, sserie = other[key], self[key]
                # index -1 does not work on series
                o_n, s_n = other.shape[1], self.shape[1]

                if key in self.residues:
                    if not chk_ceil(oserie[o_n]):
                        msg = ('Last item of {} column does not match the'
                               ' ceil {}.')
                        return FailDetail(msg.format(key, ceil))
                else:
                    if not chk_tol(sserie[s_n], oserie[o_n]):
                        msg = ('Last item of {} column does not match the'
                               ' tolerance {}.')
                        return FailDetail(msg.format(key, tol))
            if abs(s_n - o_n) > tol_iter:
                return FailDetail('Difference between number of iteration'
                                  ' is above tol_iter')
            return True

else:
    yaml_not_available_tag('Table', 'Pandas module is not available')
    yaml_not_available_tag('GwSigmaData', 'Pandas module is not available')
    yaml_not_available_tag('EtotIters', 'Pandas module is not available')
