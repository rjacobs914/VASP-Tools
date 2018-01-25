from VASP_Analyzer import OutcarAnalyzer, DirectoryUtilities
import os

oa = OutcarAnalyzer()
directory_list = DirectoryUtilities().get_downmost_directory_list()

for directory in directory_list:
    os.chdir(directory)
    print "Analyzing bandgap in directory", directory
    bandgap = oa.get_bandgap()
    f = open("bandgap.txt", "w")
    print bandgap
    f.write(str(bandgap)+"\n")
    f.close()
