#! /usr/bin/env python3
# -*- coding: utf-8 -*-


"""
    Usage:
        main.py FILE [--naccess PATH] [--points NUM] [--resolution RES] [--slice SLICE]

    Options:
        -h, --help                   Show this
        -n PATH, --naccess PATH      Absolute path to local naccess binary
        -p NUM, --points NUM         Number of points to generate on the
                                        hemisphere to criss-cross the protein.
                                        A high number will give better results,
                                        but longer calculations. [default: 250]
        -r RES, --resolution RES     An integer in angströms setting the
                                        resolution of the slicing. Setting a
                                        resolution to 1 will create slices along
                                        a processed line with a step of 1 angström.
                                        A higher resolution will create less
                                        slices of the protein, and will then
                                        reduce the computing time [default: 5].
        -s SLICE, --slice SLICE      An integer in angströms setting the thickness
                                        of a slice. A thick slice will potentially
                                        englobe more accessible residues and reduce
                                        the number of slices. This reduces the
                                        computations but is less resolutive. [default: 15]
"""


# IMPORTS

from Bio.PDB import NACCESS
from Bio.PDB import PDBParser
from docopt import docopt
from collections import abc
from operator import itemgetter
from datetime import datetime
from multiprocessing import Pool, cpu_count
from functools import partial
from mpl_toolkits.mplot3d import Axes3D
from itertools import product
from copy import deepcopy

import copy
import math
import tracemalloc
import sys

import numpy as np
import matplotlib.pyplot as plt
import src.protein as protein
import src.sphere as sphere
import src.vector as vector
import src.pdb as pdb


def loop(processed_lines, prot_dict, thickness, resolution, sphere_point):
    """Main loop doing the calculations for 1 line of the sphere.

        Args:
            processed_lines: An array of dictionaries containing the processed lines
            prot_dict: A dictionary compiling coordinates and more infos on
                        every c_alphas of the protein
            sphere_point: Coordinates of a point of the hemisphere.
                            Iterative argument of type Vector(x, y, z). Changes
                            on each iteration of the function for the parallelization

        Returns:
            The line processed during the current iteration,
            appended to the array containing all other lines processed
    """
    # print("Point de la sphère: ", sphere_point)
    # We set a plane far away from the protein (70 angströms)
    plane_normal = vector.Vector(sphere_point[0], sphere_point[1], sphere_point[2]) * 70
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
    # and a 15 angströms (by default) parallel plane.
    # The slice is sliding along the plane's normal towards the farthest
    # c_alpha with a step of 5 (by default) angströms each loop
    nb_steps = math.ceil((longest_distance - shortest_distance) / resolution)
    # Dictionary containing the current processed line's informations:
    #   - Average hydrophobicity
    #   - Relative hydrophobicity in all slices
    #   - Number of steps
    #   - Shortest distance between residues and the plane
    line = {"slice_hydro": np.zeros((nb_steps, 1)), # The index of the array = slice_step
            "line_average_hydro": None,
            "nb_steps": nb_steps,
            "shortest_distance": shortest_distance}
    # print("Number of slices to treat for the current line: ", nb_steps)
    # Process all the slices of the current line / direction
    process_slices_of_line(dist_ca_to_plane, line, nb_steps, thickness, resolution)
    # Calculate an average of the hydrophobicity of the line
    line_average_hydro = np.sum(line["slice_hydro"]) * nb_steps**2
    assert isinstance(line_average_hydro, (float, int)), "Error 5: The line average hydrophobicity should be an int or float"
    # Stock the plane's normal vector and the average hydrophobicity in list as tuples:
    # Example: [ (plane_normal_1, average_hydrophobicity_1), ... ]
    line["line_average_hydro"] = tuple((plane_normal, line_average_hydro))
    processed_lines.append(line)
    # print("Hydrophobicity of the slices in line:\n\tMin: {:.4f}\n\tMax: {:.4f}".format(np.amin(line["slice_hydro"]), np.amax(line["slice_hydro"])))
    # print("Average hydrophobicity of the line: \nPlane normal vector: ", plane_normal, "\nAverage hydrophobicity: {:.4f}".format(line_average_hydro), "\n")
    return processed_lines


def process_slices_of_line(dist_ca_to_plane, line, nb_steps, thickness, resolution):
    """Calculate the relative hydrophobicity of all slices of the current line.
        The slice slides along its line with a step of 1 angström

        Args:
            dist_ca_to_plane: Distances between all c_alphas and the plane.
                                Distances are being updated for each iteration.
            line: Dictionary compiling the line's slices hydrophobicity info and
                    the average_hydrophobicity of the line
            nb_steps: Number of slices to create along the current line to cover
                        the whole protein
        Returns:
            Nothing. Updates values inplace.
    """
    # 1 step = 1 slice
    for step in range(nb_steps):
        residues_in_slice = []
        # Check which c_alphas are in the 15 angströms slice
        for res_id, dist in dist_ca_to_plane.items():
            if 0.0 <= dist <= thickness:
                residues_in_slice.append(prot_dict[res_id]["resName"])
        # Calculate the relative hydrophobicity of the slice
        line["slice_hydro"][step] = protein.slice_relative_hydrophobicity(residues_in_slice, len(residues_in_slice))
        # Slide the plane along the line. Step = 5 angströms by default
        dist_ca_to_plane = { res_id: dist - resolution for res_id, dist in dist_ca_to_plane.items() }
    return None


if __name__ == '__main__':

    # For runtime stat
    startTime = datetime.now()

    # Parse command line
    arguments = docopt(__doc__, version='Transmembrane Protein Areas 1.0')
    pdb_file = arguments["FILE"]
    thickness = int(arguments["--slice"])
    resolution = int(arguments["--resolution"])

    # Run NACCESS with the Biopython wrapper
    pdb_struct = PDBParser(QUIET=True)
    struct = pdb_struct.get_structure(pdb_file[:4].upper(), pdb_file)
    model = struct[0]

    # Use custom naccess installation path if specified in command line argument
    if arguments["--naccess"]:
        rsa_data, asa_data = NACCESS.run_naccess(model, pdb_file,
                                                 naccess=arguments["--naccess"])
    else:
        rsa_data, asa_data = NACCESS.run_naccess(model, pdb_file)
    # Parse the naccess output .rsa file to retrieve
    # the relative % of solvant accessible area for each CA
    naccess_rsa = NACCESS.process_rsa_data(rsa_data)
    # Keep only residues having a relative accessibility > 30 (arbitrary)
    accessible_residues = protein.keep_accessible_residues(naccess_rsa)

    # The dict compiles informations on the protein residues.
    # It contains solvant accessible c_alpha coordinates,
    # the residues names and their respective solvant accessibility area value
    # We also get the center of mass of the protein.
    prot_dict, center_of_mass = protein.build_prot_dict(
        pdb_file, accessible_residues)

    # Generate n points on a hemisphere englobing the protein. By default n = 250.
    if arguments["--points"]:
        nb_points = int(arguments["--points"])
        sphere_points = sphere.generate_points_on_sphere(int(nb_points))
    else:
        nb_points = 250 # Default
        sphere_points = sphere.generate_points_on_sphere(nb_points)

    # We scale coordinates to be in a (0, 0, 0) centered
    # coordinates system to simplify further calculations
    # The coordinates are modified inplace.
    prot_dict = protein.scale_ca_coords(prot_dict, center_of_mass)


    ########################################
    # Main calculations loop is parallelized
    ########################################

    # List containing the processed sphere lines
    processed_lines = []
    # Parallelization of the main loop
    # The calculations for each line of the hemisphere is parallelized.
    # This means that several lines are calculated simultaneously
    pool = Pool(processes=cpu_count())
    func = partial(loop, processed_lines, prot_dict, thickness, resolution)
    processed_lines = pool.imap(func, sphere_points)
    pool.close()
    pool.join()
    # Extract the "best" line, the one maximizing the average hydrophobicity
    # processed_lines = protein.get_best_slice(processed_lines)
    best_results = protein.get_best_results(processed_lines)

    print("Best line: \n\t", "Point of the sphere: ", best_results[0][0] + center_of_mass, "\n\tHydrophobicity: {:.4f}".format(best_results[0][1]), "\n")
    print("Center of mass: ", center_of_mass)

    # We generate the points simulating the membranes
    pts_mb_1, pts_mb_2 = protein.generate_membranes(processed_lines, best_results, resolution)

    new_pdb_file = pdb.write_pdb(pdb_file, pts_mb_1, pts_mb_2)
    # Write a small PyMol script to visualize the best line
    pdb.write_pml_script(best_results[0][0] + center_of_mass, new_pdb_file)

    ################# RUNTIME STATS ##################################
    ##################################################################
    print("\n\nProgram runtime: ", datetime.now() - startTime)
