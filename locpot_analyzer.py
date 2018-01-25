#!/share/apps/EPD_64bit/bin/python2.7
"""
__author__ = 'Ryan Jacobs'
__version__ = '1.0'
__date__ = 'Last updated September 27, 2016'
__email__ = 'rjacobs3@wisc.edu'

This app is designed to make it easy to perform plotting and extract useful features from electrostatic potential
information stored in the LOCPOT file.

To generate the LOCPOT file, you must run VASP jobs with the following setting in the INCAR file:

LVTOT = True

If you are running a surface calculation and are interested in calculating the surface potential (work function), you
will also need the following tags. Note that it is assumed your surface termination is oriented in the c-direction, which
is considered standard practice and is the result the surface_maker app will also give you.

IDIPOL = 3
LDIPOL = True

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

# Specify whether you'd like to generate a plot of the planar-averaged (x-y plane) electrostatic potential vs z coord (True/False)
plot_electrostatic_potential = True

# Specify whether you'd like to calculate the work function. Only relevant for surface slabs (True/False)
calculate_work_function = True

# Specify whether you'd like to use the Helmholtz equation to calculate the empirically expected difference in
# work function between top and bottom surface based on the VASP calculated dipole. (True/False)
calculate_delta_work_function_helmholtz = False

########################################################################################################################
########################################################################################################################

from VASP_PostProcessing import LocpotAnalyzer
from VASP_Analyzer import OutcarAnalyzer
import logging

def main():
    # Create log file
    logging.basicConfig(filename='locpot_analysis.log', level=logging.INFO)

    loc = LocpotAnalyzer(poscar="POSCAR", outcar="OUTCAR", locpot="LOCPOT")
    oa = OutcarAnalyzer(outcar="OUTCAR")

    if plot_electrostatic_potential == bool(True):
        loc.get_electrostatic_potential()
        logging.info('The electrostatic potential energy plot and data have been saved!')

    if calculate_work_function == bool(True):
        (work_function_top, work_function_bot, vacuum_energy_top, vacuum_energy_bot) = loc.get_workfunction()
        fermi_energy = oa.get_fermi_energy()
        # Write results to log
        logging.info('The fermi energy of the system is %3.3f' % (fermi_energy))
        logging.info('The vacuum energy of the top surface is %3.3f eV' % (vacuum_energy_top))
        logging.info('The vacuum energy of the bottom surface is %3.3f eV' % (vacuum_energy_bot))
        logging.info('The work function of the top surface is %3.3f eV' % (work_function_top))
        logging.info('The work function of the bottom surface is %3.3f eV' % (work_function_bot))

    if calculate_delta_work_function_helmholtz == bool(True):
        delta_workfunction = loc.get_empirical_delta_workfunction()
        logging.info('The Helmholtz equation gives an empirical work function difference between top and bottom surface of %3.3f eV' % (delta_workfunction))

if __name__=="__main__":
    main()

