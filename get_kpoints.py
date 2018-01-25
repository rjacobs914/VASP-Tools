import os
from VASP_Setup import *
from VASP_Analyzer import DirectoryUtilities

cwd = os.getcwd()
dir_list = DirectoryUtilities().get_downmost_directory_list()
for dir_entry in dir_list:
    print "Entering directory", dir_entry
    os.chdir(dir_entry)
    kpt_setup = KpointsFileSetup(mesh_type="Monkhorst-Pack", mesh_value="2 2 2")
    kpt_setup.write_kpoints_file()
    #if "HSE" in dir_entry:
    #    kpt_setup = KpointsFileSetup(mesh_type="Monkhorst-Pack", mesh_value="1 1 1")
    #    kpt_setup.write_kpoints_file()
    #else:
    #    kpt_setup = KpointsFileSetup(mesh_type="Monkhorst-Pack", mesh_value="4 4 4")
    #    kpt_setup.write_kpoints_file()
