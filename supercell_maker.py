#!/share/apps/EPD_64bit/bin/python2.7
"""
__author__ = 'Ryan Jacobs'
__version__ = '1.0'
__date__ = 'Last updated September 27, 2016'
__email__ = 'rjacobs3@wisc.edu'

This app is designed to make it very easy to create supercells of any material contained in the Materials Project (or
from an existing POSCAR file.

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

scaling_matrix = [[2, 0, 0], [0, 2, 0], [0, 0, 2]]

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

if __name__ == "__main__":
    main()
