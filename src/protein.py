"""
.. module:: protein
  :synopsis: This module contains functions that process PDB files to retrieve
               3D coordinates of alpha carbons and calculate the center of mass of the protein.
               The parsing is made only once to optimize calculations.

.. moduleauthor:: Gabriel Cretin M2 BIB
"""

from src.vector import *
from operator import itemgetter
import sys



def get_com(x, y, z, nb_ca):
    """Calculate the Center Of Mass from a list of coordinates

    Args:
        x: Cumulative sum of C_alpha's x coordinates
        y: Cumulative sum of C_alpha's y coordinates
        z: Cumulative sum of C_alpha's z coordinates
        nb_ca: Number of alpha carbons in the protein

    Returns:
        src.vector.Vector: The center of mass of the protein as Vector(x, y, z)
    """
    return Vector(x, y, z) / nb_ca

def keep_accessible_residues(naccess_rsa):
    """From the output of naccess we keep only accessible residues
    which have a all_atoms_rel value > 30 (arbitrary threshold)

    Args:
        naccess_rsa: A dictionnary containing the output of naccess's calculations

    Returns:
        dict: Keys are the residue ids and as value their solvant accessible area
    """
    accessible_residues_dict = {}
    for (chain_id, res_id), data_dict in naccess_rsa.items():
        for key, val in data_dict.items():
            if key == "all_atoms_rel" and val >= 30:
                accessible_residues_dict[res_id[1]] = val
    return accessible_residues_dict


def build_prot_dict(pdb_file, accessible_residues):
    """1. Get the coordinates of alpha carbones in the PDB
       2. Check if the residue is accessible to solvant
       3. Builds a dictionnary compiling infos on:
            - 3D coordinates of accessible residues
            - value of relative accessibility to solvant
            - residue name

        Args:
            pdb_file: The protein's PDB file
            accessible_residues: Dictionnary containing the accessible residues and their relative accessibility value

        Returns:
            dict: Residue_id:
                    - Vector(x, y, z)
                    - All_atoms_rel_accessibility_value
                    - Residue_name

            Vector: Protein's center_of_mass
    """
    prot_dict = {}
    nb_ca = 0
    x_com = 0
    y_com = 0
    z_com = 0
    with open(pdb_file, 'r') as file_in:
        for line in file_in:
            atom_type = line[0:6].strip()  # "ATOM " or "HETATM"
            atom_name = line[12:16].strip()
            if atom_type == "ATOM" and atom_name == "CA":
                residue_name = str(line[17:20].strip())
                residue_num = int(line[22:26].strip())
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                # Cumulative sum of coordinates for the calculation
                # of the center of mass
                x_com += x
                y_com += y
                z_com += z
                nb_ca += 1
                # Keep the residue if it is accessible to solvant
                # Build a dictionnary compiling all infos for the residue
                if residue_num in accessible_residues:
                    prot_dict[residue_num] = {'3Dcoords': Vector(x, y, z),
                                              'all_atoms_rel': accessible_residues[residue_num],
                                              'resName': residue_name}
    return (prot_dict, get_com(x_com, y_com, z_com, nb_ca))


def scale_ca_coords(prot_dict, center_of_mass):
    """Place the cartesian system centered in (0, 0, 0) origin

        Args:
            prot_dict: Coordinates of all c_alphas of the protein
        """
    for c_alpha, infos in prot_dict.items():
        infos['3Dcoords'] = infos['3Dcoords'] - center_of_mass
    return prot_dict


def slice_relative_hydrophobicity(residues, nb_residues_in_slice):
    """Calculates the relative hydrophobicity of a list of residues

        Returns:
            float: :math:`Relative\ hydrophobicity = \\frac{hydrophobe\ residues}{total\ residues}`

    """
    hydrophobes = ["PHE", "ILE", "GLY", "LEU", "MET", "TRP", "TYR", "VAL"]
    try:
        rel_hydro = len(set(residues).intersection(hydrophobes)) / nb_residues_in_slice
    except ZeroDivisionError as err:
        sys.exit("It seems like there is no residues in the slice: " + str(err))
    assert isinstance(rel_hydro, (float, int)) == True, "Error 3: rel_hydro should be an int or a float."
    return rel_hydro



def max_sub_array_sum(array):
    """Implementation of the Kadane's algorithm to solve the maximum sub-array
    problem in O(n) time and O(1) space. It finds the maximum contiguous subarray
    and print its starting and end index. Here it will return the indexes of the
    slices between which there is the maximum hydrophobicity: the area of the membrane !

        Args:
            array: Numpy array containing the relative hydrophobicity of all
                    the slices of the best line

        Returns:
            tuple: (a, b, c) such that sum(array[a:b]) == c and c is maximal
    """
    best = current = 0
    current_index = start_index = best_index = 0
    for index, val in enumerate(array):
        if current + val > 0:
            current += val
        else: # reset start position
            current, current_index = 0, index + 1
        if current > best:
            start_index, best_index, best = current_index, index + 1, current
    return (start_index, best_index, best)



def get_best_results(processed_lines):
    """Parse the results of the parallelization.
    Finds the line with the highest average hydrophobicity value and determines
    the range of slices maximizing the accessible residues hydrophobicity,
    which will point to the famous transmembrane area !

        Args:
            processed_lines: A list containing lines (dictionaries) with their
                                respective slices infos and values of average
                                hydrophobicity
        Returns:
            list: [(plane_normal, average_hydrophobicity), start_index,
                    best_index, best, nb_steps, shortest_distance]
    """
    lines = []
    best_slices = None
    nb_steps = None
    shortest_distance = None
    # The result of the parallelization is an ImapIterator,
    # which explains the necessity of doing the double for loop.
    for index in processed_lines:
        for line in index:
            lines.append(line)
    # Get the best line
    best_line = max([line["line_average_hydro"] for line in lines], key=itemgetter(1))
    # Get the slices of the best line
    for line in lines:
        if line["line_average_hydro"] == best_line:
            best_slices = line["slice_hydro"]
            nb_steps = line["nb_steps"]
            shortest_distance = line["shortest_distance"]
    # Get the indexes of the slices between which there is the maximum hydrophobicity
    start_index, best_index, best = max_sub_array_sum(best_slices)
    return [best_line, start_index, best_index, best, nb_steps, shortest_distance]


def generate_membranes(processed_lines, best_results, resolution):
    """Generate points in the space to simulate the membranes.
    They will be represented in PyMol at the end as 2 planes like both membranes.

        Args:
            processed_lines: A list containing lines (dictionaries) with their
                                respective slices infos and values of average
                                hydrophobicity
            best_results: A list [(plane_normal, average_hydrophobicity),
                                    start_index, best_index, best, nb_steps,
                                    shortest_distance]
        Returns:
            tuple: points_membrane_1, points_membrane_2
    """
    # Retrieve the results
    plane_normal = best_results[0][0]
    shortest_distance = best_results[5]
    nb_steps = best_results[4]
    start_index = best_results[1]
    best_index = best_results[2]

    # Distance of membrane 1 to best plane
    dist_m1 = resolution * (start_index + 1)
    # Distance of membrane 2 to best plane
    dist_m2 = resolution * (best_index + 1)

    plane_normal = np.array([plane_normal.x, plane_normal.y, plane_normal.z])

    # We generate 2 * 500 dummy points simulating the membranes
    point1  = np.array([plane_normal[0] + shortest_distance + dist_m1, plane_normal[1] + shortest_distance + dist_m1, plane_normal[2] + shortest_distance + dist_m1])
    point2  = np.array([plane_normal[0] + shortest_distance + dist_m2, plane_normal[1] + shortest_distance + dist_m2, plane_normal[2] + shortest_distance + dist_m2])

    # a plane is a*x+b*y+c*z+d=0
    # [a,b,c] is the plane_normal. Thus, we have to calculate
    # d and we're set
    d1 = -point1.dot(plane_normal)
    d2 = -point2.dot(plane_normal)

    # create x,y
    X, Y = np.meshgrid(range(20), range(20))
    positions = np.vstack([X.ravel(), Y.ravel()])

    # calculate corresponding z
    z1 = (-plane_normal[0] * positions[1] - plane_normal[1] * positions[0] - d1) * 1. / plane_normal[2]
    z2 = (-plane_normal[0] * positions[1] - plane_normal[1] * positions[0] - d2) * 1. / plane_normal[2]

    points_membrane_1 = points_membrane_2 = np.zeros((400, 3))
    points_membrane_1[:, 0] = positions[1]
    points_membrane_1[:, 1] = positions[0]
    points_membrane_1[:, 2] = z1
    points_membrane_2[:, 0] = positions[1]
    points_membrane_2[:, 1] = positions[0]
    points_membrane_2[:, 2] = z2

    return points_membrane_1, points_membrane_2
