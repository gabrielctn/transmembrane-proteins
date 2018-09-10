

def get_COM(list_coords):
    """Calculate the Center Of Mass of list of coordinates"""
    assert isinstance(
        list_coords, list) == True, "Error 1: Input should be a list !"
    assert len(list_coords) != 0, "Error 1: The list is empty !"
    nb_coords = float(len(list_coords))
    x = 0.0
    y = 0.0
    z = 0.0
    for coord_atom in list_coords:
        x += coord_atom['x']
        y += coord_atom['y']
        z += coord_atom['z']
    return [x / nb_coords, y / nb_coords, z / nb_coords]


def get_ca_coords(pdb_file):
    """Get the coordinates of alpha carbones in the PDB"""
    list_CA = []
    with open(pdb_file, 'r') as f:
        for line in f:
            atom_type = line[0:6].strip()  # "ATOM " or "HETATM"
            atom_name = line[12:16].strip()
            chain = line[21:22].strip()
            if atom_type == "ATOM" and atom_name == "CA":
                residue_name = str(line[17:20].strip())
                residue_num = int(line[22:26].strip())
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                list_CA.append({'resNum': residue_num, 'x': x,
                                'y': y, 'z': z, 'resName': residue_name})
    return list_CA
