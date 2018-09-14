#! /usr/bin/env python3
# -*- coding: utf-8 -*-


"""
    Usage:
        main.py FILE [--naccess PATH] [--points NUM]

    Options:
        -h --help             Show this
        -n --naccess  PATH    Absolute path to local naccess binary
        -p --points   NUM     Number of points to generate on the hemisphere to criss-cross the protein.
                                A high number will give better results, but longer calculations. [default: 250]
"""


# IMPORTS

from Bio.PDB import NACCESS
from Bio.PDB import PDBParser
from docopt import docopt
from collections import abc
from operator import itemgetter

import numpy as np
import copy
import math

import src.pdb as pdb
import src.sphere as sphere
import src.vector as vector


if __name__ == '__main__':
    # Parse command line
    arguments = docopt(__doc__, version='Transmembrane Protein Areas 1.0')
    ###
    ### TODO: check_args(arguments)
    ###

    print(arguments)
    pdb_file = arguments["FILE"]

    # Run NACCESS with the Biopython wrapper
    pdb_struct = PDBParser()
    struct = pdb_struct.get_structure(pdb_file[:4].upper(), pdb_file)
    model = struct[0]
    # Use custom naccess installation path if specified in command line argument
    if arguments["--naccess"] == True:
        rsa_data, asa_data = NACCESS.run_naccess(model, pdb_file,
                                                 naccess=arguments["PATH"])
    else:
        rsa_data, asa_data = NACCESS.run_naccess(model, pdb_file)
    # Parse the naccess output .rsa file to retrieve
    # the relative % of solvant accessible area for each CA
    naccess_rsa = NACCESS.process_rsa_data(rsa_data)
    # Keep only residues having a relative accessibility > 30 (arbitrary)
    accessible_residues = pdb.keep_accessible_residues(naccess_rsa)
    # print(accessible_residues)

    # The dict compiles informations on the protein residues.
    # It contains solvant accessible c_alpha coordinates,
    # the residues names and their respective solvant accessibility area value
    # We also get the center of mass of the protein.
    prot_dict, center_of_mass = pdb.build_prot_dict(
        pdb_file, accessible_residues)

    # Generate n points on a hemisphere englobing the protein. By default n = 250.
    if arguments["--points"]:
        sphere_points = sphere.generate_points_on_sphere(arguments["NUM"])
    else:
        sphere_points = sphere.generate_points_on_sphere(250)
    # We scale coordinates to be in a (0, 0, 0) centered
    # coordinates system to simplify further calculations
    # The coordinates are modified inplace.
    prot_dict = pdb.scale_ca_coords(prot_dict, center_of_mass)
    nb_residues = len(prot_dict)

    ########################
    # Main calculations loop
    ########################
    for point in sphere_points:
        print("Point de la sphère: ", point)
        # We set a plane far away from the protein (1000 angströms)
        plane_normal = vector.Vector(point[0], point[1], point[2]) * 1000
        # Search the nearest and farthest c_alpha to the plane
        # We keep all distances between c_alphas and the plane in a dictionary
        # to simplify calculations when sliding the plane along the line.
        dist_ca_to_plane = {}
        for c_alpha, infos in prot_dict.items():
            dist_ca_to_plane[c_alpha] = infos['3Dcoords'].dist_to_plane(plane_normal)
        min_res_id, shortest_distance = min(dist_ca_to_plane.items(), key=itemgetter(1))
        max_res_id, longest_distance = max(dist_ca_to_plane.items(), key=itemgetter(1))
        # We place the plane at the level of the nearest c_alpha.
        # For this we subtract the shortest distance to all other
        # c_alpha - plane distances
        dist_ca_to_plane = {res_id: dist - shortest_distance for res_id, dist in dist_ca_to_plane.items()}
        # Calculate relative hydrophobicity of all c_alphas between this plane
        # and a 15 angströms parallel plane.
        # The slice is sliding along the plane's normal towards the farthest
        # c_alpha with a step of 1 angström each loop
        rel_hydrophobicity = 0
        nb_steps = math.ceil(longest_distance - shortest_distance)
        print("Nombre de tranches à traiter pour cette droite: ", nb_steps)
        for step in range(1, nb_steps):
            residues_in_slice = []
            # Check which c_alphas are in the 15 angströms slice
            for res_id, dist in dist_ca_to_plane.items():
                if 0.0 <= dist <= 15.0:
                    residues_in_slice.append(prot_dict[res_id]["resName"])
            # Calculate the relative hydrophobicity of the slice
            rel_hydrophobicity += pdb.relative_hydrophobicity(residues_in_slice, nb_residues)
            # The plane slides of step 1 towards the furthest c_alpha
            dist_ca_to_plane = {res_id: dist - 1 for res_id, dist in dist_ca_to_plane.items()}
        # for k, v in dist_ca_to_plane.items():
        #     print(k, v)
        average_hydrophobicity = rel_hydrophobicity / nb_steps
        print("Somme d'hydrophobicité relative:", rel_hydrophobicity)
        print("Hydrophobicité moyenne de la droite:", average_hydrophobicity)
