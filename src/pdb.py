from shutil import copyfile
import os



def generate_membranes_dum():
    pass


def write_pdb(pdb_file, plane_normal):
    """ Copy the original pdb into a new file and write the coordinates of DUM
    atoms to visualize the membranes
    """
    # Copy the original pdb file to a new pdb_file_new.pdb
    base, ext = os.path.splitext(pdb_file)
    new_file = base + "_new" + ext
    command = "cp " + pdb_file + " " + new_file
    os.system(command)
    # Write the new lines at the end of the file
    with open(new_file, "w") as file_in:
        line = "{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}".format("HETATM", "DUM" )
        file_in.write()


def write_pml_script(sphere_point):
    """Write a small PyMol script to visualize the best line"""
    with open("src/pymol_visualize.pml", "w") as file_in:
        line1 = "pseudoatom pt1, pos=[{}, {}, {}]".format(-sphere_point.x,
                                                                -sphere_point.y,
                                                                -sphere_point.z)
        line2 = "pseudoatom pt2, pos=[{}, {}, {}]".format(sphere_point.x,
                                                                sphere_point.y,
                                                                sphere_point.z)
        line3 = "distance /pt1, /pt2"
        line4 = "set dash_gap, 0"
        line5 = "set dash_radius, 0.3"
        line6 = "set dash_round_ends, 0"
        line7 = "set dash_color, 0xffcc00, dist01"
        line8 = ("hide labels, dist01")
        file_in.write('{}\n{}\n{}\n{}\n{}\n{}\n{}\n{}'.format(line1,line2,line3,
                                                              line4,line5,line6,
                                                              line7, line8))
