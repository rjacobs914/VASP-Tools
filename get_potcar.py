import os
from VASP_Setup import *
from VASP_Analyzer import DirectoryUtilities

cwd = os.getcwd()
dir_list = DirectoryUtilities().get_downmost_directory_list()
for dir_entry in dir_list:
    print "Entering directory", dir_entry
    os.chdir(dir_entry)
    ptr_setup = PotcarFileSetup(psp_type="PBE", poscar="POSCAR")
    ptr_setup.write_potcar_file()
