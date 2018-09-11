#! /usr/bin/env python3
# -*- coding: utf-8 -*-


"""Usage:
    main.py FILE

   Options:
    -h --help   Show this
"""


# IMPORTS

from Bio.PDB import NACCESS
from Bio.PDB import PDBParser
from docopt import docopt
import process_pdb as ppdb
import sphere



if __name__ == '__main__':
    # Parse command line
    arguments = docopt(__doc__, version='Transmembrane Protein Areas 1.0')
    pdb_file = arguments["FILE"]

    # Run NACCESS with Biopython wrapper
    pdb_struct = PDBParser()
    struct = pdb_struct.get_structure(pdb_file[:4].upper(), pdb_file)
    model = struct[0]
    # Set the probe_size to 1.0 angstrom to stick to paper's method.
    rsa_data, asa_data = NACCESS.run_naccess(model, pdb_file, "1")
    naccess_rel_dict = NACCESS.process_rsa_data(rsa_data)
    naccess_atom_dict = NACCESS.process_asa_data(asa_data)

    # Extract C_alpha coordinates
    list_CA = ppdb.get_ca_coords(pdb_file)
    # Calculate Centre of Mass
    center_of_mass = ppdb.get_com(list_CA)
    # print(list_CA)
    # print(center_of_mass)

    sphere_points = generate_points_on_sphere(center_of_mass, num_points)
    
