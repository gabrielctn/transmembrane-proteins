"""
.. module:: protein
  :synopsis: This module contains functions that process PDB files to retrieve
               3D coordinates of alpha carbons and calculate the center of mass of the protein.
               The parsing is made only once to optimize calculations.

.. moduleauthor:: Gabriel Cretin M2 BIB
"""

from src.vector import *
from operator import itemgetter



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
            if key == "all_atoms_rel" and val >= 20:
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
    rel_hydro = len(set(residues).intersection(hydrophobes)) / nb_residues_in_slice
    assert isinstance(rel_hydro, (float, int)) == True, "Error 3: rel_hydro should be an int or a float."
    return rel_hydro


def get_best_line(processed_lines):
    """Finds the line with the highest average hydrophobicity value

        Args:
            processed_lines: A list containing lines (dictionaries) with their
                                respective normal vector and values of average
                                hydrophobicity
        Returns:
            tuple: (plane_normal, average_hydrophobicity)
    """
    tuples = []
    # The result of parallelization returns a nested list,
    # which explains the double for loop
    for index in processed_lines:
        for line in index:
            tuples.append(line["line_average_hydro"])
    return max(tuples, key=itemgetter(1))


def get_best_slice(processed_lines):
    """This function finds the slice maximizing the hydrophobicity of the
    accessible residues inside"""
    best_slices = None
    for index in processed_lines:
        for line in index:
            print("coucou")
    return processed_lines
    #         if best_line == line["line_average_hydro"]:
    #             best_slices = line["slice_hydro"]
    # # Get the index of the maximal hydrophobicity value among all slices
    # index_max = np.argmax(best_slices)
    # print("Best slices: \n", best_slices)
    # print("Index, value of best slice: ", index_max, best_slices[index_max])




def maxSubArraySum(array, size):
    """Function to find the maximum contiguous subarray
    and print its starting and end index"""
    max_so_far = -maxsize - 1
    max_ending_here = 0
    start = 0
    end = 0
    s = 0

    for i in range(0,size):
        max_ending_here += array[i]

        if max_so_far < max_ending_here:
            max_so_far = max_ending_here
            start = s
            end = i

        if max_ending_here < 0:
            max_ending_here = 0
            s = i+1

    print ("Maximum contiguous sum is %d"%(max_so_far))
    print ("Starting Index %d"%(start))
    print ("Ending Index %d"%(end))




def max_sub_array_sum(list):
    """Implementation of the Kadane's algorithm to solve the maximum sub-array problem
    in O(n) time and O(1) space

    Returns:
        tuple: (a, b, c) such that sum(list[a:b]) == c and c is maximal"""
    best = current = 0
    current_index = start_index = best_index = 0
    for index, val in enumerate(array):
        if current + val > 0:
            current += val
        else: # reset start position
            current, current_index = 0, index + 1
        if current > best:
            start_index, best_index, best = current_index, index + 1, current
    return start_index, best_index, best




# best = cur = 0
#     for i in l:
#         cur = max(cur + i, 0)
#         best = max(best, cur)
#     return best
