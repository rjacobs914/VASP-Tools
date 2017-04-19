#!/share/apps/EPD_64bit/bin/python2.7
"""
__author__ = 'Ryan Jacobs'
__version__ = '1.0'
__date__ = 'Last updated September 27, 2016'
__email__ = 'rjacobs3@wisc.edu'

This app is designed to make it easy to perform plotting and extract useful features from the density of states
information contained in the DOSCAR file.

To generate the correct DOS information, you must run your VASP jobs with the following tags in your INCAR file:

ISPIN = 2
LORBIT = 11

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

# Specify whether or not to create a plot of the Total Density of States (True/False)
plot_total_dos = True

# Specify whether or not to create a plot of the Projected Density of States (True/False)
plot_projected_dos = True

# Specify whether or not to calculate the bandgap (True/False)
calculate_bandgap = True

# If calculate_bandgap = True, need to specify whether to get bandgap from DOSCAR or OUTCAR (or both)
get_bandgap_from_dos = False
get_bandgap_from_outcar = True

# Specify whether or not to calculate the band centers of each element (True/False)
calculate_bandcenters = True

# Specify whether or not to calculate the oxygen charge transfer gap (only relevant for systems containing oxygen) (True/False)
calculate_charge_transfer_gap = True

# Specify whether or not to calculate the effective densities of states of the valence and conduction bands
calculate_effective_dos = True

########################################################################################################################
########################################################################################################################

from VASP_PostProcessing import *

dos = DOSAnalyzer(poscar="POSCAR", incar="INCAR", outcar="OUTCAR", doscar="DOSCAR")

def main():
    if plot_total_dos == True:
        dos.get_total_dos_plot
    if plot_projected_dos == True:
        dos.get_projected_dos_plot
    if calculate_bandcenters == True:
        dos.get_bandcenters
    if calculate_bandgap == True:
        if get_bandgap_from_dos == True:
            dos.get_bandgap_from_dos
        if get_bandgap_from_outcar == True:
            dos.get_bandgap_from_outcar
    if calculate_charge_transfer_gap == True:
        dos.get_chgtransgap
    if calculate_effective_dos == True:
        dos.get_effective_dos

if __name__=="__main__":
    main()
