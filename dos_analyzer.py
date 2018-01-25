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
plot_total_dos = False

# Specify whether or not to create a plot of the Projected Density of States (True/False)
plot_projected_dos = True
include_plot_labels = False
pickle_plot = False
use_custom_colors = True

# Specify whether or not to calculate the bandgap (True/False)
calculate_bandgap = True

# Specify whether or not to calculate the band centers of each element (True/False)
calculate_bandcenters = False
write_bandcenters_to_file = False

# Specify whether or not to calculate the oxygen charge transfer gap (only relevant for systems containing oxygen) (True/False)
calculate_charge_transfer_gap = False

# Specify whether or not to calculate the effective densities of states of the valence and conduction bands
calculate_effective_dos = False
effective_dos_temperature = 1000 # Temperature, in K, to evaluate effective DOS at

########################################################################################################################
########################################################################################################################

from VASP_PostProcessing import DoscarAnalyzer
from matplotlib import pyplot
dos = DoscarAnalyzer(poscar="POSCAR", incar="INCAR", outcar="OUTCAR", doscar="DOSCAR")

def main():
    if plot_total_dos == bool(True):
        dos.plot_total_dos()
    if plot_projected_dos == bool(True):
        plot = dos.plot_projected_dos(plot_labels=str(include_plot_labels), pickle_plot=str(pickle_plot), use_custom_color_type=str(use_custom_colors))
    if calculate_bandcenters == bool(True):
        dos.get_bandcenters(write_dicts_to_file=write_bandcenters_to_file)
    if calculate_bandgap == bool(True):
        dos.get_bandgap_from_dos()
    if calculate_charge_transfer_gap == bool(True):
        dos.get_O_chargetransfergap()
    if calculate_effective_dos == bool(True):
        dos.get_effective_dos(temperature=effective_dos_temperature)

if __name__=="__main__":
    main()
