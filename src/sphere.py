"""
.. module:: sphere
  :synopsis: This module contains the function that generates of points
                evenly distributed at the surface of a hemisphere.

.. moduleauthor:: Gabriel Cretin M2 BIB
"""

import numpy as np
from src.vector import *


def generate_points_on_sphere(num_points):
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

    # Translate the spherical coordinates to cartesian 3D coordinates (x, y, z)
    # Set the radius to 100 angstroms to englobe most of proteins
    points = np.zeros((num_points, 3))
    points[:, 0] = np.cos(theta) * np.sin(phi)
    points[:, 1] = np.sin(theta) * np.sin(phi)
    points[:, 2] = np.cos(phi)
    return points
