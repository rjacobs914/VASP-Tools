import os
from VASP_Setup import *
from VASP_Analyzer import DirectoryUtilities

cwd = os.getcwd()
dir_list = DirectoryUtilities().get_downmost_directory_list()
for dir_entry in dir_list:
    print "Entering directory", dir_entry
    os.chdir(dir_entry)
    inc_setup = IncarFileSetup()
    inc_setup.write_predefined_incar_file(simulation_type="full_relax", xc_functional="SCAN", number_of_nodes=2, disable_symmetry=False, poscar="POSCAR",
                                    material_type="insulator", write_chgcar=False, write_wavecar=False)
