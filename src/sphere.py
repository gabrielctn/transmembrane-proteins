import numpy as np
from src.vector import *

"""
.. module:: sphere
   :synopsis: This module contains the function that generates of points
   evenly distributed at the surface of a hemisphere.

.. moduleauthor:: Gabriel Cretin M2 BIB
"""


def generate_points_on_sphere(com_coordinates, num_points):
    """Generate *num_points* points evenly distributed on a hemisphere,
    centered on the center of mass of the protein thanks to the
    golden angle 3 - sqrt(5).

    Args:
        com_coordinates: Coordinates of the center of mass of the protein
        num_points: Number of desired points on the hemisphere

    Returns:
        A Numpy array of n 3D cartesian coordinates 
        centered on the center of mass of the protein
    """
    # Generate n indices on an evenly distributed interval
    indices = np.arange(0, num_points, dtype=float) + 0.5
    golden_angle = np.pi * (3 - np.sqrt(5))
    # We need only points on an hemisphere instead of whole sphere.
    phi = np.arccos(1 - indices / num_points)
    theta = golden_angle * indices

    # Translate the spherical coordinates to cartesian 3D (x, y, z)
    # Center the hemisphere to the center of mass of the protein
    # The radius of the hemisphere is 1.
    points = np.zeros((num_points, 3))
    points[:, 0] = com_coordinates.x + np.cos(theta) * np.sin(phi)
    points[:, 1] = com_coordinates.y + np.sin(theta) * np.sin(phi)
    points[:, 2] = com_coordinates.z + np.cos(phi)
    return points
