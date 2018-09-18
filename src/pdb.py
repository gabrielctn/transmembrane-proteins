"""This module implements functions to visualize the program's results
with PyMol. It generates a PDB file and a .pml file."""


from shutil import copyfile
import os


def write_pdb(pdb_file, points_membrane_1, points_membrane_2):
    """ Copy the original pdb into a new file and write the coordinates of DUM
    atoms to visualize the membranes

        Args:
            pdb_file: Path to the pdb_file
            points_membrane_1: 3D coordinates of the points representing the first membrane
            points_membrane_2: 3D coordinates of the points representing the second membrane

        Returns:
            file: New PDB file with DUM atoms simulating the membranes
    """
    # Copy the original pdb file to a new pdb_file_new.pdb
    base, ext = os.path.splitext(pdb_file)
    new_pdb_file = base + "_new" + ext
    command = "cp " + pdb_file + " " + new_pdb_file
    os.system(command)
    # Write the new lines at the end of the file
    with open(new_pdb_file, "a") as file_in:
        i = 10000
        for point in points_membrane_1:
            if i % 2 == 0:
                line = "{:6s}{:5d}      {:3s} {:1s} {:5d}       {:8.3f} {:8.3f} {:8.3f}\n".format(
                    "HETATM ", i, "N", "DUM", i, float(point[0]), float(point[1]), float(point[2]))
            else:
                line = "{:6s}{:5d}      {:3s} {:1s} {:5d}       {:8.3f} {:8.3f} {:8.3f}\n".format(
                    "HETATM ", i, "O", "DUM", i, float(point[0]), float(point[1]), float(point[2]))
            i += 1
            file_in.write(line)
        i = 10000
        for point in points_membrane_2:
            if i % 2 == 0:
                line = "{:6s}{:5d}      {:3s} {:1s} {:5d}       {:8.3f} {:8.3f} {:8.3f}\n".format(
                    "HETATM ", i, "N", "DUM", i, float(point[0]), float(point[1]), float(point[2]))
            else:
                line = "{:6s}{:5d}      {:3s} {:1s} {:5d}       {:8.3f} {:8.3f} {:8.3f}\n".format(
                    "HETATM ", i, "O", "DUM", i, float(point[0]), float(point[1]), float(point[2]))
            i += 1
            file_in.write(line)
    return new_pdb_file


def write_pml_script(sphere_point, pdb_file):
    """Write a small PyMol script to visualize the best line/direction

        Args:
            sphere_point: A point on the hemisphere englobing the protein
            pdb_file: Path to the pdb_file
    """
    with open("src/pymol_visualize.pml", "w") as file_in:
        line1 = "cmd.load('{}')".format(pdb_file)
        line2 = "cmd.pseudoatom('pt1', pos=[{}, {}, {}])".format(-sphere_point.x,
                                                                 -sphere_point.y,
                                                                 -sphere_point.z)
        line3 = "cmd.pseudoatom('pt2', pos=[{}, {}, {}])".format(sphere_point.x,
                                                                 sphere_point.y,
                                                                 sphere_point.z)
        line4 = "cmd.distance('/pt1', '/pt2')"
        line5 = "cmd.set('dash_gap', '0')"
        line6 = "cmd.set('dash_radius', '0.3')"
        line7 = "cmd.set('dash_round_ends', '0')"
        line8 = "cmd.set('dash_color', '0xffcc00', 'dist01')"
        line9 = "cmd.hide('labels', 'dist01')"
        file_in.write('{}\n{}\n{}\n{}\n{}\n{}\n{}\n{}\n{}'.format(line1, line2, line3,
                                                                  line4, line5, line6,
                                                                  line7, line8, line9))
