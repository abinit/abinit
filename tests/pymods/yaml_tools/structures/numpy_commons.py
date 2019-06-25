'''
    Define structures depending on Numpy to be found in Abinit YAML formated
    output extending the possible operations on the extracted data.
'''
from __future__ import print_function, division, unicode_literals
import numpy as np
from ..common import BaseArray
from ..register_tag import yaml_seq

yaml_seq(BaseArray)


@yaml_seq
class Atoms3D(BaseArray):
    pass


@yaml_seq
class CartForces(Atoms3D):
    pass


@yaml_seq
class Matrix33(BaseArray):
    '''
        Define a matrix of shape (3, 3) compatible with numpy arrays and
        with YAML tags.
    '''
    def __init__(self, shape=(3, 3), *args, **kwargs):
        assert shape == (3, 3)
        super(Matrix33, self).__init__(shape, *args, **kwargs)

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
