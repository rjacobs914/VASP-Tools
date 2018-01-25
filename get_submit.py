import os
from VASP_Setup import *
from VASP_Analyzer import DirectoryUtilities

cwd = os.getcwd()
dir_list = DirectoryUtilities().get_downmost_directory_list()
for dir_entry in dir_list:
    print "Entering directory", dir_entry
    os.chdir(dir_entry)
    sfs = SubmitFileSetup(cluster="ACI", queue="morgan2", number_of_nodes=2)
    sfs.write_submit_file()
