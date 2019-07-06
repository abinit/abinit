'''
    Define structures depending on numpy to be found in Abinit YAML formatted
    output extending the possible operations on the extracted data.
'''
from __future__ import print_function, division, unicode_literals
from ..common import BaseArray
from ..register_tag import yaml_seq

yaml_seq(BaseArray)


@yaml_seq
class Atoms3D(BaseArray):
    '''Base class for (natom, 3) arrays.'''


@yaml_seq
class CartForces(Atoms3D):
    '''Cartesian forces as (natom, 3) array'''


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
class TensorCart(Matrix33):

    def is_symetric(self):
        '''Return true if tensor is symmetric.'''
        for i in range(3):
            for j in range(i, 3):
                if self[i, j] != self[j, i]:
                    return False
        return True

    def is_anti_symetric(self):
        '''Return true if tensor is anti-symmetric '''
        for i in range(3):
            for j in range(i, 3):
                if self[i, j] != -self[j, i]:
                    return False
        return True
