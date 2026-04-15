"""
Define structures depending on numpy to be found in Abinit YAML formatted
output extending the possible operations on the extracted data.
"""
from ..common import BaseArray
from ..register_tag import yaml_seq

yaml_seq(BaseArray)


@yaml_seq
class Atoms3D(BaseArray):
    """Base class for (natom, 3) arrays."""


@yaml_seq
class CartForces(Atoms3D):
    """Cartesian forces as (natom, 3) array"""


@yaml_seq
class Matrix33(BaseArray):
    """Define a matrix of shape (3, 3) compatible with numpy arrays and with YAML tags."""

    def __init__(self, shape=(3, 3), *args, **kwargs):
        """
        Initialize the 3x3 matrix.

        Args:
            shape (tuple, optional): Shape of the matrix. Must be (3, 3).
            *args: Variable length argument list.
            **kwargs: Arbitrary keyword arguments.
        """
        assert shape == (3, 3)
        super().__init__(shape, *args, **kwargs)

    @classmethod
    def from_seq(cls, s):
        """
        Create a Matrix33 from a sequence.

        Args:
            s (sequence): The sequence containing matrix elements.

        Returns:
            Matrix33: A new 3x3 matrix.
        """
        new = super().from_seq(s)
        assert new.shape == (3, 3)
        return new


@yaml_seq
class CartTensor(Matrix33):

    def is_symmetric(self, tol_abs=1e-8):
        """
        Check if the tensor is symmetric.

        Args:
            tol_abs (float, optional): Absolute tolerance. Defaults to 1e-8.

        Returns:
            bool: True if symmetric, False otherwise.
        """
        for i in range(3):
            for j in range(i, 3):
                if abs(self[i, j] - self[j, i]) > tol_abs: return False
        return True

    def is_antisymmetric(self, tol_abs=1e-8):
        """
        Check if the tensor is anti-symmetric.

        Args:
            tol_abs (float, optional): Absolute tolerance. Defaults to 1e-8.

        Returns:
            bool: True if anti-symmetric, False otherwise.
        """
        for i in range(3):
            for j in range(i, 3):
                if abs(self[i, j] + self[j, i]) > tol_abs: return False
        return True
