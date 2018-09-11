"""This module contains all the functions used for the search 
of the membrane pane. That includes the generation of points 
evenly distributed at the surface of a sphere, and lines 
passing by the center of mass and these points.
"""

def generate_points_on_sphere(com_coordinates, num_points):
    """Generate *num_points* points evenly distributed on a sphere centered 
    on the center of mass of the protein using the Vogel's method
    with the golden angle.
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
    points[:,0] = com_coordinates[0] + radius * np.cos(theta)
    points[:,1] = com_coordinates[1] + radius * np.sin(theta)
    points[:,2] = com_coordinates[2] + fixed_z_axis
    return points


def get_vector(com_coordinates, sphere_point):
    """Returns the vector between the 
    center of mass to the point on the sphere
    """
    return np.substract(sphere_point, com_coordinates)


def get_vector_magnitude(vect):
    """Calculate the magnitude (Euclidian norm) of a vector"""
    return np.sqrt(vect.dot(vect)) # Faster than np.linalg.norm(x)


def get_unit_vector(vect):
    """Normalize a vector by its magnitude: return the unit vector"""
    return np.true_divide(vect, get_vector_magnitude(vect))


def distance_point_from_plane(point, plane):
    return (p.A*v.x + p.B * v.y + p.C* v.z + p.D) / np.sqrt(p.A*p.A + p.B*p.B + p.C*p.C) 