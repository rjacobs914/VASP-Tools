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

# If you want to analyze the stability of a material you simulated with DFT, set to True. If set to False, need to
# specify your own material composition and energy
analyze_data_from_VASP_run = True

# If analyze_data_from_VASP_run = False, specify values for the arguments below:
custom_composition_dictionary = None # e.g., [{"La":2, "O":3} , {"La":1, "Mn":1, "O":3}]
custom_energy_list = None # e.g., [-100, -90, -526]

# If you want to remove certain materials from the Materials Project database, so they aren't considered in the stability
# analysis, provide a list of their material id's here

material_ids_to_remove = None # e.g. ["mp-123", "mp-12345"]

# If you just want to use DFT energies and not incorporate other specific chemical potential values, set this to False.
# If you want to specify your own chemical potentials of gaseous species to correspond to some type of environmental
# conditions, set this as True and also set the below values to specify how the chemical potential is calculated

use_my_own_chemical_potentials = True

# If you have specified use_my_own_chemical_potentials = True above, you must set the additional following parameters

# Set which species you want to calculate the chemical potential of, in the form of a list of strings. For example,
# ["O", "H"]

species_of_chemical_potential = ["O", "H"]

# Set the temperature (in K) of the environment to calculate the chemical potential at

temperature = 1073

# Set the partial pressure (in atm) of the species to calculate the chemical potential at

partial_pressure_of_species = 0.2

########################################################################################################################
# If you are calculating the chemical potential of H specifically, there are a couple other variables to set:

# Specifies whether or not to calculate H chemical potential as being in equilibrium with H2O. True = equilibrium with
# H2O. False = H2 gas

calc_H_chem_pot_H2Oequil = True

# If specifying calc_H_chem_pot_H2Oequil = True, you must specify what relative humidity (RH) of H2O vapor you want.
# 0 = 0% RH, 1 = 100% RH.

relative_humidity_H2Oequil = 1

# If specifying calc_H_chem_pot_H2Oequil = True, you must specify the pressure of O2 gas (in atm) which is also in equilibrium
# with the water. Typically, this will be the partial pressure of O2 in air, which is 0.2 atm

oxygen_pressure_H2Oequil = 0.2

########################################################################################################################
# If performing calculations where you want a set of organic molecule energies considered, set this to true
include_organic_molecules = False

# If include_organic_molecules = True, decide if you want to apply energy shift to molecules
include_organic_molecule_energy_shift = False

########################################################################################################################
########################################################################################################################

from VASP_PostProcessing import StabilityAnalyzer, ChemicalPotentialAnalyzer
import os

#mapi_key = os.environ['MAPI_KEY']
mapi_key = "TtAHFCrZhQa7cwEy"

def main():
    if use_my_own_chemical_potentials == bool(True):
        custom_chem_pot_dict = {}
        for entry in species_of_chemical_potential:
            if entry == "O":
                O_chem_pot, H, S = ChemicalPotentialAnalyzer(temperature=temperature, pressure=partial_pressure_of_species,
                                                       functional="PBE", energy_shift=False).get_O_chem_pot()
                custom_chem_pot_dict["O"]=O_chem_pot
            if entry == "H":
                if calc_H_chem_pot_H2Oequil == True:
                    H_chem_pot = ChemicalPotentialAnalyzer(temperature=temperature, pressure=partial_pressure_of_species,
                                                           functional="PBE", energy_shift=False, relative_humidity=relative_humidity_H2Oequil,
                                                           pressure_O2_forH2O=oxygen_pressure_H2Oequil).get_H_chem_pot_fromH2O()
                    custom_chem_pot_dict["H"]=H_chem_pot
                if calc_H_chem_pot_H2Oequil == False:
                    H_chem_pot = ChemicalPotentialAnalyzer(temperature=temperature, pressure=partial_pressure_of_species,
                                                           functional="PBE", energy_shift=False).get_H_chem_pot()
                    custom_chem_pot_dict["H"]=H_chem_pot

    elif use_my_own_chemical_potentials == bool(False):
        custom_chem_pot_dict = None

    stability = StabilityAnalyzer(mapi_key=mapi_key, poscar="POSCAR", oszicar="OSZICAR", get_data_from_VASP_files=analyze_data_from_VASP_run,
                 additional_elements_to_include=species_of_chemical_potential, material_ids_to_remove=material_ids_to_remove,
                                      composition_dict=custom_composition_dictionary, composition_energy=custom_energy_list)
    stability.get_phase_diagram(use_custom_chem_pots=use_my_own_chemical_potentials, custom_chem_pot_dict=custom_chem_pot_dict,
                                    include_organic_molecules=include_organic_molecules, include_organic_molecule_shift=include_organic_molecule_energy_shift)

if __name__=="__main__":
    main()