import numpy as np
from numbers import Number

"""
.. module:: vector
   :synopsis: This module implements a Vector class and functions associated to
   Vectors (coordinates) manipulation

.. moduleauthor:: Gabriel Cretin M2 BIB
"""


class Vector:
    """
    .. class:: Vector
    This class implements 3D vectors supporting numpy arrays

    Attributes:
        x: x coordinate
        y: y coordinate
        z: z coordinate
        coords: Numpy array containing the 3D coordinates
    """

    def __init__(self, x=0, y=0, z=0):
        """ Creates a vector, example: v = Vector(1, 2, 3)
        with individual coordinates + a numpy array
        """
        self.x = float(x)
        self.y = float(y)
        self.z = float(z)
        self.coords = np.array([x, y, z])

    def norm(self):
        """Calculate the norm (length, magnitude) of a vector"""
        return np.sqrt(self.coords.dot(self.coords))  # Faster than np.linalg.norm()

    def unit(self):
        """Normalize a vector by its magnitude: returns the unit vector"""
        return np.true_divide(self.coords, self.coords.norm)

    def vector_to(self, point):
        """Returns the vector between self and the point given in argument"""
        if isinstance(point, Vector):
            return point.coords - self.coords
        return point - self.coords

    def distance_to_plane(self, normal):
        """Calculates the distance between a 3D point and a plane

        Args:
            normal: the normal vector of the plane
        """
        # a plane is a*x+b*y+c*z+d=0
        # [a,b,c] is the normal. Thus, we have to calculate d
        d = -(normal.norm()**2)
        return np.true_divide(np.abs(normal.x * self.x +
                                     normal.y * self.y +
                                     normal.z * self.z +
                                     d), np.sqrt(normal.x ** 2 +
                                                 normal.y ** 2 +
                                                 normal.z ** 2))

    def create_line_between(self, dest_point):
        pass


    def create_far_plane(com, norm, ):
        pass



    def __str__(self):
        """String representation"""
        return "<{}, {}, {}>".format(self.x, self.y, self.z)

    def __neg__(self):
        """Signing"""
        return Vector(-self.x, -self.y, -self.z)

    def __mul__(self, scalar):
        """Scalar Multiplication"""
        return Vector(self.x * scalar,
                      self.y * scalar,
                      self.z * scalar)

    def __div__(self, scalar):
        """Division"""
        return np.true_divide(self.__copy().coords, scalar)

    def __add__(self, operand):
        """Addition"""
        return Vector(self.x + operand.x,
                      self.y + operand.y,
                      self.z + operand.z)

    def __sub__(self, operand):
        """Subtraction"""
        return self.__copy() + -operand

    def __getitem__(self, index):
        if isinstance(index, int):
            return self.coords[index]
        # slicing
        elif isinstance(index, type(slice(1))):
            return self.coords.__getitem__(index)
        else:
            raise KeyError

    def __setitem__(self, key, value):
        if not isinstance(value, numbers.Number):
            raise ValueError
        if (key == 0):
            self.x = value
        elif (key == 1):
            self.y = value
        elif (key == 2):
            self.z = value
