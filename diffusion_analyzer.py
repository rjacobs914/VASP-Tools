#!/share/apps/EPD_64bit/bin/python2.7
"""
__author__ = 'Ryan Jacobs'
__version__ = '1.0'
__date__ = 'Last updated September 27, 2016'
__email__ = 'rjacobs3@wisc.edu'

This app is designed to make it easy to compute and plot the mean-squared displacement from Ab-initio Molecular
Dynamics (AIMD) runs (data in the XDATCAR file) and also calculate the resulting diffusion coefficient of a given species.
You can combine multiple XDATCAR files together (cat XDATCAR1 XDATCAR2 > XDATCAR) to conduct the analysis over a longer
simulated time. In the future, there will be some option added to have the code do this for you automatically.

To conduct an AIMD run, you must run VASP with the following settings in your INCAR file:

ALGO = Fast
PREC = Low
POTIM = 2 (time step in femtoseconds, 2 is a reasonable choice)
SMASS = 0
NSW = 10000 (number of time steps, set based on total desired time and time step size)
IBRION = 0
TEBEG = 1000 (beginning temperature, have same as TEEND if doing constant T)
TEEND = 1000 (ending temperature, have same as TEBEG if doing constant T)
NBLOCK = 1

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

# Specify which element specie in your calculation you'd like to calculate the diffusion coefficient of (e.g. "N")
element_of_diffusing_species = "O"

# Specify the size of the timestep (in femtoseconds) used in your AIMD run. This should match the POTIM tag in the INCAR.
timestep_size = 2

# Specify the temperature (in K) used in the AIMD run.
temperature = 2073

# Specify which method will be used to smooth out the MSD. Must choose either "max" or "constant". "max" is recommended.
msd_smoothing_method = "max"

# Specify whether you'd like to create and save a plot of the mean squared displacement. Useful to analyze to see if
# diffusion is really occurring.
plot_mean_squared_displacement = True

# Specify names of the XDATCAR file(s) you wish to analyze. Files from multiple AIMD runs should be placed in the same
# directory with names as e.g. XDATCAR1, XDATCAR2, etc.
xdatcar_file_names = ['XDATCAR']

########################################################################################################################
########################################################################################################################

from VASP_PostProcessing import DiffusionAnalyzerAIMD

def main():
    diff = DiffusionAnalyzerAIMD(xdatcar=xdatcar_file_names, diffusing_species=element_of_diffusing_species, temperature=temperature,
                        timestep=timestep_size, steps_to_ignore=100, smoothing=msd_smoothing_method, min_obs=30,
                        avg_nsteps=1000, plot_msd=plot_mean_squared_displacement)
    diff.get_diffusion_analysis_from_xdatcars()

if __name__=="__main__":
    main()
