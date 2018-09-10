

def generate_points_on_sphere(com_coordinates):
    """ http://blog.marmakoide.org/?p=1 """
    num_points = 256 # Number of points we want
    golden_angle = np.pi * (3 - np.sqrt(5))
    theta = golden_angle * np.arange(n)
    z = np.linspace(1 - 1.0 / num_points, 1.0 / num_points - 1, num_points)
    radius = np.sqrt(1 - z * z)

    # Generate points centered on the center of mass
    com_x, com_y, com_z = com_coordinates
    points = np.full((n, 3), [com_x, com_y, com_z])
    points[:,0] = radius * np.cos(theta)
    points[:,1] = radius * np.sin(theta)
    points[:,2] = z
    return points
