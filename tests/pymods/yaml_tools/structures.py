'''
    Define some structures to be found in Abinit YAML formated output
    extending the possible operations on the extracted data.
'''
from __future__ import print_function, division, unicode_literals
from .abinit_metadata import ITERATOR_RANKS
from .register_tag import (yaml_map, yaml_seq, yaml_auto_map,
                           yaml_implicit_scalar, yaml_scalar)
import numpy as np


@yaml_map
class IterStart(object):
    '''
        Mark the begining of a iteration of a given iterator.
    '''
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


@yaml_seq
class AutoNumpy(np.ndarray):
    '''
        Define a base class for YAML tags converted to numpy compatible
        objects.  Can be used for converting any YAML array of number of any
        dimension into a numpy compatible array.
    '''
    __yaml_tag = 'Array'

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
class CartForces(AutoNumpy):
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
        new = AutoNumpy.from_seq(s)
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


@yaml_scalar
class WhiteSpaceSeparated(object):
    __yaml_tag = 'WSS'

    def __init__(self, column_labels, data):
        self.columns = column_labels
        self.data = data

    @classmethod
    def from_scalar(cls, scal):
        header, lines = scal.split('\n')
        data = []
        for line in lines:
            data.append([])
            for v in line.split():
                if v.isdigit():
                    data[-1].append(int(v))
                else:
                    try:
                        data[-1].append(float(v))
                    except ValueError:
                        data[-1].append(v)
        return cls(header, data)

    def to_scalar(self):
        length = max(1 + max(self.header, key=lambda x: len(x)), 8)
        line_fmt = ("{:>" + str(length) + "}") * len(self.header)
        scal = line_fmt.format(*self.header) + '\n'
        for line in self.data:
            scal += line_fmt.format(*line) + '\n'
        return scal

    # TODO create a practical interface compatible with python 2 and 3
    # def __getitem__(self, key):
    #     if isinstance(key, int):
    #         return WhiteSpaceSeparated(self.columns,
    #                                    self.data[slice(key, key + 1)])
    #     if isinstance(key, str):
    #         j = self.columns.index(str)
    #         return [
    #             line[j] for line in self.data
    #         ]
    #     if isinstance(key, slice):
    #         if isinstance(key.start, int):
    #             return WhiteSpaceSeparated(self.columns, self.data[key])
    #         elif isinstance(key.start, str):
    #             j = self.columns.index(key.start)
    #             k = self.columns.index(key.stop)
    #             hview = self.columns[j:k]
    #             view = [
    #                 line[j:k] for line in self.data
    #             ]
    #             return WhiteSpaceSeparated(hview, view)

    def __iter__(self):
        for i, column in self.columns:
            yield column, [
                line[i] for line in self.data
            ]


@yaml_auto_map
class AbinitError(object):
    __yaml_tag = 'ERROR'
    pass


@yaml_auto_map
class AbinitWarning(object):
    __yaml_tag = 'WARNING'
    pass


@yaml_auto_map
class AbinitInfo(object):
    __yaml_tag = 'INFO'
    pass


@yaml_auto_map
class AbinitComment(object):
    __yaml_tag = 'COMMENT'
    pass
