#!/share/apps/EPD_64bit/bin/python2.7
"""
__author__ = 'Ryan Jacobs'
__version__ = '1.0'
__date__ = 'Last updated September 27, 2016'
__email__ = 'rjacobs3@wisc.edu'

This app is designed to make it very easy to create surfaces of any material contained in the Materials Project
database using simple user input.

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

# Specify whether you'd like to use an existing POSCAR to generate a surface or import the structure from Materials Project

use_existing_poscar = True

# If you specified "False" above, you must specify a material composition below (e.g. 'CaTiO3').
# Specifying "True" above overrides this entry
material = 'W'

# The scaling matrix to create a supercell from the conventional cell of your material. To use just the conventional
# cell, use the identity matrix as the scaling matrix.

scaling_matrix = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]

# Choose the surface termination you want as a set of Miller indices. The surface will be oriented such that the chosen
# miller indices forms the surface termination and will point in the c-direction

surface_miller_index = [0, 0, 1]

# Choose the minimum thickness of your surface slab, in units of Angstroms. You may need to experiment with this setting
# to obtain the right slab thickness

surface_slab_thickness = 20

# Choose the thickness of the vacuum region above the slab, in units of Angstroms. You want the vacuum to be thick
# enough so the calculation is converged with respect to properties like surface energy, work function, etc. Typically
# 20 is a good value

vacuum_region_thickness = 20

# Choose whether or not to shift the atom coordinates of the slab to the bottom of the supercell. Setting to True will
# shift atoms to bottom of supercell, which can sometimes help with surface run convergence if dipole corrections
# are used

shift_slab_to_bottom = True

#################################################################################################
# The following entries are optional if you want to remove layers from your surface slab
#################################################################################################

# Choose whether you want to remove any layers from the surface slab you've created. Can be useful for isolating certain
# surface terminations or making your slab the correct desired size.
trim_layers_from_slab = False

# Specify the number of top layers to remove
number_of_top_layers_to_remove = 0

# Specify the number of bottom layers to remove
number_of_bottom_layers_to_remove = 0

#################################################################################################
# The following entries are optional if you want to apply Selective Dynamics to your surface slab
#################################################################################################

# Set this to True if you want to use the Selective Dynamics functionality. Otherwise a POSCAR will be created with
# no Selective Dynamics

apply_selective_dynamics = True

# This setting dictates whether Selective Dynamics will be applied to both the top and bottom surfaces. If set to False,
# Selective Dynamics will only be applied to the top surface.

relax_top_and_bottom = True

# This setting determines the number of material layers from the surface to be relaxed. It is up to the user to know
# based on their own situation and testing what a reasonable value for this parameter is.

number_of_layers_to_relax = 2

########################################################################################################################
########################################################################################################################

from VASP_Setup import *
import os
import warnings

def main():
    mapi_key = os.environ['MAPI_KEY']
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        if use_existing_poscar == bool(False):
            make_psr = PoscarFileSetup(material_composition=material, mapi_key=mapi_key, get_only_most_stable_structure=True)
            make_psr.write_poscar_file()
    mod_psr = PoscarFileModifier(poscar="POSCAR")
    mod_psr.make_poscar_supercell(scaling_matrix=scaling_matrix)
    mod_psr.make_poscar_slab(surface_hkl=surface_miller_index, minimum_thickness=surface_slab_thickness,
                             minimum_vacuum=vacuum_region_thickness, shift_slab_to_bottom=shift_slab_to_bottom)
    if trim_layers_from_slab == bool(True):
        surface_trimming_dict = {"Top": number_of_top_layers_to_remove, "Bottom": number_of_bottom_layers_to_remove}
        mod_psr.trim_poscar_slab(surface_trimming_dict=surface_trimming_dict)
    if apply_selective_dynamics == bool(True):
        mod_psr.add_selective_dynamics(number_of_layers=number_of_layers_to_relax, top_and_bottom=relax_top_and_bottom)

if __name__ == "__main__":
    main()
