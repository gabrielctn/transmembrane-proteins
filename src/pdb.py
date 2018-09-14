from src.vector import *


"""
.. module:: pdb
   :synopsis: This module contains functions that process PDB files to retrieve
   3D coordinates of alpha carbons and calculate the center of mass of the protein.
   The parsing is made only once to optimize calculations.

.. moduleauthor:: Gabriel Cretin M2 BIB
"""


def get_com(x, y, z, nb_ca):
    """Calculate the Center Of Mass from a list of coordinates

    Args:
        x: Cumulative sum of C_alpha's x coordinates
        y: Cumulative sum of C_alpha's y coordinates
        z: Cumulative sum of C_alpha's z coordinates
        nb_ca: Number of alpha carbons in the protein

    Returns:
        tuple: The center of mass of the protein as (x, y, z)

    """
    return Vector(x / nb_ca, y / nb_ca, z / nb_ca)


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
            - dict: prot_dict -> residue_id -> Vector(x, y, z)
                                          -> all_atoms_rel_accessibility_value
                                          -> residue_name
            - Vector: Protein's center_of_mass

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
    """Place the cartesian system centered in (0, 0, 0) origin"""
    for c_alpha, infos in prot_dict.items():
        infos['3Dcoords'] = infos['3Dcoords'] - center_of_mass
    return prot_dict
