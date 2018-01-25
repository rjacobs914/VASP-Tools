from VASP_Analyzer import DirectoryUtilities, PoscarAnalyzer
import os
import shutil

directory_list = DirectoryUtilities().get_downmost_directory_list()

for directory in directory_list:
    os.chdir(directory)
    print "In directory", directory
    incar = open("INCAR", "r")
    incarnew = open("INCARnew", "w")

    for line in incar:
        if ("NELMIN" not in line):
            incarnew.write(line)
    #incarnew.write("ALGO = GW0"+"\n")
    #incarnew.write("NELM = 1"+"\n")
    #incarnew.write("NPAR = 2"+"\n")

    #if "BSCF2.625" in directory:
    #    incarnew.write("NUPDOWN = 28"+"\n")
    #if "GBCO" in directory:
    #    incarnew.write("NUPDOWN = 48"+"\n")
    #if "LaCoO3" in directory:
    #    incarnew.write("NUPDOWN = 16"+"\n")
    #if "LaCrO3" in directory:
    #    incarnew.write("NUPDOWN = 24"+"\n")
    #if "LaCuO3" in directory:
    #    incarnew.write("NUPDOWN = 8"+"\n")
    #if "LaFeO3" in directory:
    #    incarnew.write("NUPDOWN = 40"+"\n")
    #if "LaMnO3" in directory:
    #    incarnew.write("NUPDOWN = 32"+"\n")
    #if "LaNiO3" in directory:
    #    incarnew.write("NUPDOWN = 8"+"\n")
    #if "LaScO3" in directory:
    #    incarnew.write("NUPDOWN = 0"+"\n")
    #if "LaTiO3" in directory:
    #    incarnew.write("NUPDOWN = 8"+"\n")
    #if "LaVO3" in directory:
    #    incarnew.write("NUPDOWN = 16"+"\n")
    #if "LSC25" in directory:
    #    incarnew.write("NUPDOWN = 18"+"\n")
    #if "LSC50" in directory:
    #    incarnew.write("NUPDOWN = 20"+"\n")
    #if "LSM25" in directory:
    #    incarnew.write("NUPDOWN = 30"+"\n")
    #if "PBCO" in directory:
    #    incarnew.write("NUPDOWN = 20"+"\n")
    #if "SBCO" in directory:
    #    incarnew.write("NUPDOWN = 20"+"\n")
    #if "SrCoO3" in directory:
    #    incarnew.write("NUPDOWN = 12"+"\n")
    #if "SrFeO2.75" in directory:
    #    incarnew.write("NUPDOWN = 36"+"\n")
    #if "SrTiO3" in directory:
    #    incarnew.write("NUPDOWN = 0"+"\n")
    #if "SrVO3" in directory:
    #    incarnew.write("NUPDOWN = 8"+"\n")
   
    #element_list = PoscarAnalyzer().get_element_names()
    #atom_amounts = PoscarAnalyzer().get_atom_amounts()

    #magmom_str = ''
    #tm_list = ["Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu"]
    #for element, atom_amount in zip(element_list, atom_amounts):
    #    if element == "Ti":
    #        magmom_str += " "+str(atom_amount)+"*"+"1"
    #    if element == "V":
    #        magmom_str += " "+str(atom_amount)+"*"+"2"
    #    if element == "Cr":
    #        magmom_str += " "+str(atom_amount)+"*"+"3"
    #    if element == "Mn":
    #        magmom_str += " "+str(atom_amount)+"*"+"4"
    #    if element == "Fe":
    #        magmom_str += " "+str(atom_amount)+"*"+"5"
    #    if element == "Co":
    #        magmom_str += " "+str(atom_amount)+"*"+"3"
    #    if element == "Ni":
    #        magmom_str += " "+str(atom_amount)+"*"+"1"
    #    if element == "Cu":
    #        magmom_str += " "+str(atom_amount)+"*"+"2"
    #    if element not in tm_list:
    #        magmom_str += " "+str(atom_amount)+"*"+"0"

    #incarnew.write("MAGMOM = "+magmom_str+"\n")

    #if "BSCF" in directory:
    #    print(element_list)
    #    print(atom_amounts)
    #    print(magmom_str)

    """
    if "HSE" in directory:
        for line in incar:
            if ("ISIF" not in line and "LWAVE" not in line):
                #print line.strip()
                incarnew.write(line)
        incarnew.write("ISIF=2"+"\n")
        incarnew.write("LWAVE=True"+"\n")
    else:
        for line in incar:
            if "ISIF" not in line:
                #print line.strip()
                incarnew.write(line)
        incarnew.write("ISIF=2"+"\n")
    """

    incarnew.write("NELMIN=13"+"\n")

    incarnew.close()
    incar.close()
    shutil.move("INCARnew", "INCAR")

