from VASP_Analyzer import VASPdata, DirectoryUtilities

directory_list = DirectoryUtilities().get_downmost_directory_list()
#custom_file_list = ["composition.txt", "O_pband_values.txt", "energy_above_hull.txt", "formation_energy.txt", "O_element_magnetization.txt"]
#custom_file_list = ["bandgap_fromdos.txt", "composition.txt", "O_pband_values.txt"]
file_list_to_string = ["composition.txt"]
custom_file_list = ["composition.txt", "O_pband_allstates.txt", "O_pband_occstates.txt", "O_pband_unoccstates.txt"]
tmlist = ["Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu"]
for tm in tmlist:
    custom_file_list.append(tm+"_dband_allstates.txt")
    custom_file_list.append(tm+"_dband_occstates.txt")
    custom_file_list.append(tm+"_dband_unoccstates.txt")
print "The custom file list is:", custom_file_list
vd = VASPdata()
vd.write_vaspdata_to_spreadsheet(directory_list=directory_list, use_custom_file_list=True, custom_file_list=custom_file_list, file_list_to_string=file_list_to_string)
