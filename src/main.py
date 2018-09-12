#! /usr/bin/env python3
# -*- coding: utf-8 -*-


"""Usage:
    main.py FILE [--naccess PATH]

   Options:
    -h --help             Show this
    -n --naccess  PATH    Absolute path to local naccess binary
"""


# IMPORTS

from Bio.PDB import NACCESS
from Bio.PDB import PDBParser
from docopt import docopt
from collections import abc
import pdb
import sphere
import vector





if __name__ == '__main__':
    # Parse command line
    arguments = docopt(__doc__, version='Transmembrane Protein Areas 1.0')
    #print(arguments)
    pdb_file = arguments["FILE"]

    # Run NACCESS with the Biopython wrapper
    pdb_struct = PDBParser()
    struct = pdb_struct.get_structure(pdb_file[:4].upper(), pdb_file)
    model = struct[0]
    # Use custom naccess installation path if specified in command line argument
    if arguments["--naccess"] == True:
        rsa_data, asa_data = NACCESS.run_naccess(model, pdb_file, naccess=arguments["PATH"])
    else:
        rsa_data, asa_data = NACCESS.run_naccess(model, pdb_file)
    # Parse the naccess output .rsa file to retrieve 
    # the relative % of solvant accessible area for each CA
    naccess_rsa = NACCESS.process_rsa_data(rsa_data)
    # Keep only residues having a relative accessibility > 30 (arbitrary)
    accessible_residues = pdb.keep_accessible_residues(naccess_rsa)
    ### print(accessible_residues)

    # Compile informations on the protein residues
    prot_dict = pdb.build_prot_dict(pdb_file, accessible_residues)
    # Centre of Mass
    center_of_mass = prot_dict['com']
    ### print(list_CA_coords)
    ### print(center_of_mass)

    # Generate 500 points on a hemisphere englobing the protein
    sphere_points = sphere.generate_points_on_sphere(center_of_mass, 500)
    #print(sphere_points)
    # # Search for nearest point to plane
    #for point in sphere_points:
        # Create a line passing by the center 
        # of mass and a point on the hemisphere

