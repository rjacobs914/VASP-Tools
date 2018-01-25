from VASP_Analyzer import PoscarAnalyzer, DirectoryUtilities
import os

pa = PoscarAnalyzer()
directory_list = DirectoryUtilities().get_downmost_directory_list()

for directory in directory_list:
    os.chdir(directory)
    print "Analyzing composition in directory", directory
    composition = pa.get_composition_dict()
    composition_str = ''
    for k, v in composition.items():
        composition_str += str(k)+str(v) 
    f = open("composition.txt", "w")
    print composition_str
    f.write(str(composition_str)+"\n")
    f.close()
