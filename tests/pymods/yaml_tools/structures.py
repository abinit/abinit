'''
    Define some structures to be found in Abinit YAML formated output
    extending the possible operations on the extracted data.
'''
from __future__ import print_function, division, unicode_literals
from .abinit_metadata import ITERATOR_RANKS
from .register_tag import is_available

if is_available:
    from .register_tag import yaml_map, yaml_seq, yaml_auto_map, yaml_implicit_scalar, yaml_scalar
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
    class DataDoc(object):
        def __init__(self, label='nothing', comment='no comment'):
            self.label = label
            self.comment = comment

    @yaml_auto_map
    class Etot(object):
        Etot_yaml_tag = 'ETOT'
        def __init__(self, label='nothing', comment='no comment'):
            self.label = label
            self.comment = comment

    @yaml_seq
    class AutoNumpy(np.ndarray):
        '''
            Define a base class for YAML tags converted to numpy compatible objects.
            Can be used for converting any YAML array of number of any dimension into
            a numpy compatible array.
        '''

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
        pass

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
            seq[0][0] = float(self[0, 0])  # we have to explicitly convert values to float
            seq[0][1] = float(self[1, 1])  # for yaml to recognise it properly
            seq[0][2] = float(self[2, 2])  # because those values are float compatible

            seq[1][0] = float(self[1, 2])  # but of a numpy custom type
            seq[1][1] = float(self[0, 2])
            seq[1][2] = float(self[0, 1])
            return seq

    @yaml_implicit_scalar
    class SelfTolerenced(float):
        yaml_pattern = r'[+-]?(\d+(\.\d*)?|\.\d+)([eE][+-]?\d+)?~(\d+(\.\d*)?|\.\d+)([eE][+-]?\d+)?'

        def __new__(cls, value, tol):
            return super().__new__(cls, value)

        def __init__(self, value, tol):
            self.tol = float(tol)

        @classmethod
        def from_scalar(cls, scal):
            val, tol = scal.split('~')
            return cls(val, tol)

        def to_scalar(self):
            return str(self) + '~' + str(self.tol)

        def __add__(self, other):
            val = self.real + other.real
            tol = self.tol
            if hasattr(other, 'tol'):
                tol += other.tol
            return SelfTolerenced(val, tol)

        def __radd__(self, other):
            val = self.real + other.real
            tol = self.tol
            if hasattr(other, 'tol'):
                tol += other.tol
            return SelfTolerenced(val, tol)

    @yaml_scalar
    class WhiteSpaceSeparated(object):
        WhiteSpaceSeparated_yaml_tag = 'WSS'

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

        def __getitem__(self, key):
            if isinstance(key, int):
                return WhiteSpaceSeparated(self.columns, self.data[slice(key, key + 1)])
            if isinstance(key, str):
                j = self.columns.index(str)
                return [
                    line[j] for line in self.data
                ]
            if isinstance(key, slice):
                if isinstance(key.start, int):
                    return WhiteSpaceSeparated(self.columns, self.data[key])
                elif isinstance(key.start, str):
                    j = self.columns.index(key.start)
                    k = self.columns.index(key.stop)
                    hview = self.columns[j:k]
                    view = [
                        line[j:k] for line in self.data
                    ]
                    return WhiteSpaceSeparated(hview, view)

        def __iter__(self):
            for i, column in self.columns:
                yield column, [
                    line[i] for line in self.data
                ]
