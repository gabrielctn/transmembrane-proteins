"""
.. module:: vector
  :synopsis: This module implements a Vector class and functions associated to
                Vectors (coordinates) manipulation

.. moduleauthor:: Gabriel Cretin M2 BIB

"""

import numpy as np
from numbers import Number



class Vector:
    """
    .. class:: Vector
      This class implements 3D vectors (x, y, z) on top of numpy arrays

    Attributes:
        x: x coordinate
        y: y coordinate
        z: z coordinate
    """

    @property
    def x(self):
        """Gives the first element of the vector"""
        return self._data[0]

    @property
    def y(self):
        """Gives the second element of the vector"""
        return self._data[1]

    @property
    def z(self):
        """Gives the third element of the vector"""
        return self._data[2]

    @x.setter
    def x(self, value):
        """Sets the first element of the vector"""
        self._data[0] = value

    @y.setter
    def y(self, value):
        """Sets the second element of the vector"""
        self._data[1] = value

    @z.setter
    def z(self, value):
        """Sets the third element of the vector"""
        self._data[2] = value

    def __init__(self, x=0, y=0, z=0):
        """ Creates a vector as numpy array
        Example: v = Vector(1, 2, 3)
        """
        self._data = np.array([x, y, z])

    def __str__(self):
        return "[{:.4f}, {:.4f}, {:.4f}]".format(self.x, self.y, self.z)

    def __repr__(self):
        return "[{:.4f}, {:.4f}, {:.4f}]".format(self.x, self.y, self.z)

    def __mul__(self, value):
        """Multiplication"""
        if isinstance(value, Vector):
            result = self._data * value._data
        elif isinstance(value, Number):
            result = self._data * value
        return Vector(result[0], result[1], result[2])

    def __rmul__(self, value):
        """Commutative version of Multiplication"""
        return self.__mul__(value)

    def __truediv__(self, value):
        """Division"""
        if isinstance(value, Vector):
            result = self._data / value._data
        elif isinstance(value, Number):
            result = self._data / value
        return Vector(result[0], result[1], result[2])

    def __rtruediv__(self, value):
        """Commutative version of Division"""
        return self.__div__(value)

    def __add__(self, value):
        """Addition"""
        if isinstance(value, Vector):
            result = self._data + value._data
        elif isinstance(value, Number):
            result = self._data + value
        return Vector(result[0], result[1], result[2])

    def __radd__(self, value):
        """Commutative version of Addition"""
        return self.__add__(value)

    def __sub__(self, value):
        """Subtraction"""
        if isinstance(value, Vector):
            result = self._data - value._data
        elif isinstance(value, Number):
            result = self._data - value
        return Vector(result[0], result[1], result[2])

    def __rsub__(self, value):
        """Commutative version of Subtraction"""
        if isinstance(value, Vector):
            result = value._data - self._data
        elif isinstance(value, Number):
            result = value - self._data
        return Vector(result[0], result[1], result[2])

    def __neg__(self):
        """Signing"""
        return Vector(-self.x, -self.y, -self.z)

    # Math functions

    def norm(self):
        """Calculate the norm (length, magnitude) of a vector"""
        norm = np.sqrt(self._data.dot(self._data)) # Faster than np.linalg.norm()
        assert np.isscalar(norm) == True, "Error 4: the norm should be a scalar."
        return norm

    def dist_to_plane(self, normal):
        """Calculates the distance between a 3D point and a plane

        Args:
            normal: The normal vector of the plane
        """
        # A plane equation is: a*x+b*y+c*z+d = 0
        # [a,b,c] is the normal. Thus, we have to calculate d
        d = -(normal.norm()**2)
        numerator = np.abs(normal.x * self.x + normal.y * self.y + normal.z * self.z + d)
        denominator = np.sqrt(normal.x ** 2 + normal.y ** 2 + normal.z ** 2)
        assert np.isscalar(numerator) == True, "Error 3: the numerator should be a scalar."
        assert np.isscalar(denominator) == True, "Error 3: the denominator should be a scalar."
        return np.true_divide(numerator, denominator)
