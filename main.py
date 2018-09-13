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


############################## IMPORTS

from Bio.PDB import NACCESS
from Bio.PDB import PDBParser
from docopt import docopt
from collections import abc
import src.pdb as pdb
import src.sphere as sphere
import src.vector as vector





if __name__ == '__main__':
    # Parse command line
    arguments = docopt(__doc__, version='Transmembrane Protein Areas 1.0')
    ###
    ###TODO: check_args(arguments)
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
    ### print(accessible_residues)

    # The dict compiles informations on the protein residues.
    # It contains solvant accessible c_alpha coordinates,
    # the residues names and their respective solvant accessibility area value
    # We also get the center of mass of the protein.
    prot_dict, center_of_mass = pdb.build_prot_dict(pdb_file, accessible_residues)

    # Generate n points on a hemisphere englobing the protein. By default 250.
    if arguments["--points"]:
        sphere_points = sphere.generate_points_on_sphere(center_of_mass, arguments["NUM"])
    else:
        sphere_points = sphere.generate_points_on_sphere(center_of_mass, 250)

    for point in sphere_points:
        # We set a plane far away from the protein
        plane_normal = 1000 * center_of_mass.vector_to(point)
        plane_normal = vector.Vector(plane_normal[0], plane_normal[1], plane_normal[2])
        # Search the nearest point (c_alpha) to the plane
        for c_alpha, infos in prot_dict.items():
            ca_plane_distance = infos['3Dcoords'].distance_to_plane(plane_normal)
            print("Normale du plan: ", plane_normal)
            print("Distance: ", ca_plane_distance)
