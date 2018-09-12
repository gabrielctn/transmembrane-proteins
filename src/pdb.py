"""This module contains functions that process PDB files to retrieve
3D coordinates of alpha carbons and calculate the center of mass of the protein.
The parsing is made only once to optimize calculations.
"""


def get_com(x, y, z, nb_ca):
    """Calculate the Center Of Mass from a list of coordinates"""
    return (x / nb_ca, y / nb_ca, z / nb_ca)


def keep_accessible_residues(naccess_rsa):
    """From the output of naccess we keep only accessible residues
    which have a all_atoms_rel value > 30 (arbitrary threshold)"""
    accessible_residues_dict = {}
    for (chain_id, res_id), data_dict in naccess_rsa.items():
        for key, val in data_dict.items():
            if key == "all_atoms_rel" and val >= 30:
                accessible_residues_dict[res_id[1]] = val
    return accessible_residues_dict


def build_prot_dict(pdb_file, accessible_residues):
    """1) Get the coordinates of alpha carbones in the PDB
       2) Check if the residue is accessible to solvant
       3) Builds a dictionnary compiling infos on:
            - Center of mass of the protein (tuple)
            - 3D coordinates of accessible residues
            - value of relative accessibility to solvant
            - residue name
        Returns:
            prot_dict -> center_of_mass
                      -> residue_id -> x
                                    -> y
                                    -> z
                                    -> all_atoms_rel_accessibility_value
                                    -> residue_name
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
                    prot_dict[residue_num] = {'x': x, 'y': y, 'z': z,
                                              'all_atoms_rel': accessible_residues[residue_num],
                                              'resName': residue_name}
    prot_dict['com'] = get_com(x, y, z, nb_ca)
    return prot_dict
