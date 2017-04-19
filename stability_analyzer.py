#!/share/apps/EPD_64bit/bin/python2.7
"""
__author__ = 'Ryan Jacobs'
__version__ = '1.0'
__date__ = 'Last updated September 27, 2016'
__email__ = 'rjacobs3@wisc.edu'

This app is designed to make it easy to perform thermodynamic phase stability analysis using the tools contained in
pymatgen. It enables the user to quickly check if their calculated material is stable using just DFT energies, or stable
under specific environmental conditions of their choosing by setting the chemical potentials of desired gaseous species
as of v1.0, only O and H are supported, but Cl, I, Br, F, and N will soon be added).

Note that it may be reasonable to consider a material "stable" if it is within about 40 meV/formula unit above the
convex hull, however a value closer to 0 is more indicative of good stability.

To obtain DFT energies that are consistent with pymatgen and thus the Materials Project, the user must use GGA-PBE type
pseudopotentials. Also, for compounds containing V, Cr, Mn, Fe, Co, Ni, or Mo, the user must use GGA+U with the following
U values (J=0 for all): {"V": "3.25", "Cr": "3.7", "Mn": "3.9", "Fe": "5.3", "Co": "3.32", "Ni": "6.2", "Mo": "4.38"}.
Fortunately, if you create your input files using the input_file_maker app and specify you're doing a GGA+U calculation,
this will all be taken care of for you.

################################################
General info to use this suite of analysis tools
################################################

You will need the following Python modules installed: VASP_Analyzer, VASP_Setup, VASP_PostProcessing, pymatgen, and
MPinterfaces

VASP_Analyzer, VASP_Setup and VASP_PostProcessing can be obtained from emailing Ryan Jacobs at rjacobs3@wisc.edu
Pymatgen can be downloaded and installed by running the following in your terminal:
    # pip install pymatgen --user
MPinterfaces can be downloaded and installed by running the following in your terminal (note that installation of
version 1.1.2 seems to fail, so version 1.2.0 is used):
    # pip install mpinterfaces=="1.2.0" --user

You will also need to set your Materials API key (obtained from Materials Project) in your .bashrc file as MAPI_KEY.
This can be done by doing the following steps at your terminal:
    # vi ~/.bashrc
    # export MAPI_KEY="abcdefg"

Lastly, you will also need to set the location of your VASP pseudopotentials in your .bashrc file as VASP_PSP_DIR.
If you have installed MAST, this should already be done for you. If not, this can be done by doing the following steps
at your terminal:
    # vi ~/.bashrc
    # export VASP_PSP_DIR="/home/user/vasp_pps"

"""

########################################################################################################################
########################################################################################################################

####################################################
# Required entries are below, with examples provided
####################################################

# If you just want to use DFT energies and not incorporate other specific chemical potential values, set this to False.
# If you want to specify your own chemical potentials of gaseous species to correspond to some type of environmental
# conditions, set this as True and also set the below values to specify how the chemical potential is calculated

use_my_own_chemical_potentials = True

# Specify the number of atoms per formula unit of your material. (E.g. for SrTiO3, number_of_atoms_per_formula_unit = 5,
# for UO2, number_of_atoms_per_formula_unit = 3)

number_of_atoms_per_formula_unit = 5

# Specify whether you want to include multiple newly calculated materials in your stability analysis. This is used to,
# for example, test whether a material is stable in the presence of competing phases. If you just want to perform stability
# analysis in the current working directory on one material, set this to False.

analyze_stability_for_multiple_materials = False

# If you want to test stability for multiple new materials at once (analyze_stability_for_multiple_materials = True above),
# then provide list of directories where materials data is located. Provide the full path in each case, e.g.
# list_of_directories_to_analyze = ["/home/user/materials_to_analyze/material1", "/home/user/materials_to_analyze/material2"]

list_of_directories_to_analyze = []

##########################################################

# If you have specified use_my_own_chemical_potentials = True above, you must set the additional following parameters

# Set which species you want to calculate the chemical potential of, in the form of a list of strings. For example,
# ["O", "H"]

species_of_chemical_potential = ["O"]

# Set the temperature (in K) of the environment to calculate the chemical potential at

temperature = 1200

# Set the partial pressure (in atm) of the species to calculate the chemical potential at

partial_pressure_of_species = 10**-9


##########################################################

# If you are calculating the chemical potential of H specifically, there are a couple other variables to set:

# Specifies whether or not to calculate H chemical potential as being in equilibrium with H2O. True = equilibrium with
# H2O. False = H2 gas

calc_H_chem_pot_H2Oequil = False

# If specifying calc_H_chem_pot_H2Oequil = True, you must specify what relative humidity (RH) of H2O vapor you want.
# 0 = 0% RH, 1 = 100% RH.

relative_humidity_H2Oequil = 1

# If specifying calc_H_chem_pot_H2Oequil = True, you must specify the pressure of O2 gas (in atm) which is also in equilibrium
# with the water. Typically, this will be the partial pressure of O2 in air, which is 0.2 atm

oxygen_pressure_H2Oequil = 0.2

########################################################################################################################
########################################################################################################################

from VASP_PostProcessing import *
import os

mapi_key = os.environ['MAPI_KEY']

def main():
    if use_my_own_chemical_potentials == True:
        custom_chem_pot_dict = {}
        for entry in species_of_chemical_potential:
            if entry == "O":
                O_chem_pot = ChemicalPotentialAnalyzer(temperature=temperature, pressure=partial_pressure_of_species,
                                                       functional="PBE", energy_shift=False).get_O_chem_pot
                custom_chem_pot_dict["O"]=O_chem_pot
            if entry == "H":
                if calc_H_chem_pot_H2Oequil == True:
                    H_chem_pot = ChemicalPotentialAnalyzer(temperature=temperature, pressure=partial_pressure_of_species,
                                                           functional="PBE", energy_shift=False, relative_humidity=relative_humidity_H2Oequil,
                                                           pressure_O2_forH2O=oxygen_pressure_H2Oequil).get_H_chem_pot_fromH2O
                    custom_chem_pot_dict["H"]=H_chem_pot
                if calc_H_chem_pot_H2Oequil == False:
                    H_chem_pot = ChemicalPotentialAnalyzer(temperature=temperature, pressure=partial_pressure_of_species,
                                                           functional="PBE", energy_shift=False).get_H_chem_pot
                    custom_chem_pot_dict["H"]=H_chem_pot

        if "O" in species_of_chemical_potential:
            include_O = True
        if "H" in species_of_chemical_potential:
            include_H = True
        stability = StabilityAnalyzer(mapi_key=mapi_key, poscar="POSCAR", oszicar="OSZICAR",
                                      atoms_per_formula_unit=number_of_atoms_per_formula_unit, always_include_O=include_O, always_include_H=include_H)
        stability.get_phase_diagram(use_custom_chem_pots=True, custom_chem_pot_dict=custom_chem_pot_dict,
                                    include_multiple_material_directories=analyze_stability_for_multiple_materials,
                                    material_directory_list=list_of_directories_to_analyze, include_organic_molecules=False,
                                    include_organic_molecule_shift=False)

    if use_my_own_chemical_potentials == False:
        stability = StabilityAnalyzer(mapi_key=mapi_key, poscar="POSCAR", oszicar="OSZICAR",
                                      atoms_per_formula_unit=number_of_atoms_per_formula_unit)
        stability.get_phase_diagram(use_custom_chem_pots=False, custom_chem_pot_dict={},
                                    include_multiple_material_directories=analyze_stability_for_multiple_materials,
                                    material_directory_list=list_of_directories_to_analyze, include_organic_molecules=False,
                                    include_organic_molecule_shift=False)

if __name__=="__main__":
    main()