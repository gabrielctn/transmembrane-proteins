#! /usr/bin/env python3
# -*- coding: utf-8 -*-


"""Usage:
    tm_prot.py FILE

   Options:
    -h --help   Show this
"""


# IMPORTS

from process_pdb import *
from docopt import docopt

###############


if __name__ == '__main__':
    # Parse command line, based on the program's usage. Cf module docstring.
    arguments = docopt(__doc__, version='Transmembrane Protein Areas 1.0')
    # print(arguments)
    list_CA = get_ca_coords(arguments["FILE"])
    center_of_mass = get_COM(list_CA)
    print(list_CA)
    print(center_of_mass)
