import numpy as np
from vector import *

"""
.. module:: sphere
   :synopsis: This module contains all the functions used for the search
   of the membrane pane. That includes the generation of points
   evenly distributed at the surface of a sphere, and lines
   passing by the center of mass and these points.

.. moduleauthor:: Gabriel Cretin M2 BIB
"""


def generate_points_on_sphere(com_coordinates, num_points):
    """Generate *num_points* points evenly distributed on a sphere centered
    on the center of mass of the protein using the Vogel's method
    with the golden angle.

    Args:
        com_coordinates: Coordinates of the center of mass of the protein
        num_points: Description of arg2

    Returns:
        A Numpy array of n points centered on the center of mass of the protein
    """
    golden_angle = np.pi * (3 - np.sqrt(5))
    theta = golden_angle * np.arange(num_points)
    # Return evenly spaced numbers over specified interval.
    # Keeping points on one hemisphere is sufficient
    # to draw lines afterwards
    fixed_z_axis = np.linspace(1 - 1.0 / num_points, 0, num_points)
    radius = np.sqrt(1 - fixed_z_axis * fixed_z_axis)

    # Generate points centered on the center of mass
    points = np.zeros((num_points, 3))
    points[:, 0] = com_coordinates[0] + radius * np.cos(theta)
    points[:, 1] = com_coordinates[1] + radius * np.sin(theta)
    points[:, 2] = com_coordinates[2] + fixed_z_axis
    return points
