#!/share/apps/EPD_64bit/bin/python2.7
"""
__author__ = 'Ryan Jacobs'
__version__ = '1.0'
__date__ = 'Last updated September 27, 2016'
__email__ = 'rjacobs3@wisc.edu'

This app is designed to generate all the needed input files to run VASP jobs. By executing this script from your
current working directory, it will loop through all sub directories, creating (or using an existing) POSCAR file and
generating the POTCAR, KPOINTS, INCAR and submit files for you, based on a number of available parameters explained
below.

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

#############
# Variables related to setting up POTCAR file
#############

# Specify the type of pseudopotentials to use. Current choices are "PBE" and "PW91"
type_of_pseudopotentials = "PBE"

#############
# Variables related to setting up KPOINTS file
#############

# Specify type of kpoint mesh. Choices are "Monkhorst-Pack", "Gamma" and "Auto". "Auto" is good if you are iterating
# through many directories, each containing different materials.
kpoint_mesh_type = "Monkhorst-Pack"

# Specify kpoint mesh. For kpoint_mesh_type = "Gamma" and "Monkhorst-Pack", this is e.g. "4 4 4", For "Auto", 25 is
# typically a good choice
kpoint_mesh_value = "4 4 4"

#############
# Variables related to setting up INCAR file
#############

# Specify the style of DFT calculation you're doing. Choices are:
    # "full_relax": Full volume + ions relax
    # "ions_relax": Fixed volume relax, but relaxes ions
    # "static_run": Single iteration, fixed volume and ions
    # "surface_run": Surface calculation with dipole corrections, local potential data
DFT_run_type = "full_relax"

# Specify the number of nodes to use in your calculation
number_of_nodes = 2

# Specify the type of exchange-correlation functional to use. Choices are:
    # "GGA": Standard GGA
    # "GGA+U": GGA with U values for transition metals consistent with Materials Project usage
    # "HSE": Hybrid GGA with standard 25% Hartree-Fock mixing
type_of_exchange_correlation_functional = "GGA"

# Specify the type of electronic structure the material has. Choices are "metal" and "insulator"
electronic_structure_type = "insulator"

# Specify whether or not to disable symmetry analysis (ISYM)
disable_symmetry = False

# Specify whether or not to write the CHGCAR file (True/False)
write_chgcar_file = False

# Specify whether or not to write the WAVECAR file (True/False)
write_wavecar_file = False

#############
# Variables related to setting up submit file
#############

# Specify the name of the cluster you're submitting the job to. Choices are "Turnbull", "ACI", "Cori", "Stampede"
cluster_name = "ACI"

# Specify the name of the queue you're submitting the job to.
    # "morgan.q" if cluster_name = "Turnbull"
    # "univ", "univ2", "morgan", or "morgan2" if cluster_name = "ACI"
    # "normal" if cluster_name = "Stampede"
    # "regular" if cluster_name = "Cori"
queue_name = "morgan"

# Whether or not to use the gamma-point VASP executable (True/False). Consider setting to True if your kpoint mesh is
# Gamma and you're using 1 kpoint (can get about 2x speed increase)
use_gamma_point_executable = False

#############
# Specify whether to submit your new jobs
#############

submit_new_jobs = False

########################################################################################################################
########################################################################################################################
import os
from VASP_Setup import *
from VASP_Analyzer import DirectoryUtilities

def main():
    mapi_key = os.environ['MAPI_KEY']
    cwd = os.getcwd()
    dir_list = DirectoryUtilities().get_downmost_directory_list()
    for dir_entry in dir_list:
        print "Entering directory", dir_entry
        os.chdir(dir_entry)
        new_dir = dir_entry + "/" + "new_run"

        if not os.path.exists(new_dir):
            os.mkdir(new_dir)

        if use_existing_poscar == bool(True):
            # Copy existing CONTCAR to new directory and make it POSCAR
            if os.path.exists(dir_entry+"/"+"CONTCAR"):
                if os.path.getsize(dir_entry+"/"+"CONTCAR") > 0:
                    print "Copying existing CONTCAR to %s" % new_dir
                    shutil.copy(dir_entry+"/"+"CONTCAR", new_dir)
                    os.chdir(new_dir)
                    shutil.move(new_dir+"/"+"CONTCAR", new_dir+"/"+"POSCAR")
                else:
                    print "No CONTCAR available. Copying existing POSCAR to %s" % new_dir
                    shutil.copy(dir_entry+"/"+"POSCAR", new_dir)
                    os.chdir(new_dir)
            else:
                print "Copying existing POSCAR to %s" % new_dir
                shutil.copy(dir_entry+"/"+"POSCAR", new_dir)
                os.chdir(new_dir)

        if use_existing_poscar == bool(False):
            psr_setup = PoscarFileSetup(material_composition=material, mapi_key=mapi_key)
            psr_setup.write_poscar_file()

        ptr_setup = PotcarFileSetup(psp_type=type_of_pseudopotentials, poscar="POSCAR")
        kpt_setup = KpointsFileSetup(mesh_type=kpoint_mesh_type, mesh_value=kpoint_mesh_value)
        inc_setup = IncarFileSetup()
        sub_setup = SubmitFileSetup(cluster=cluster_name, queue=queue_name, number_of_nodes=number_of_nodes,
                                        gamma_point=use_gamma_point_executable)
        ptr_setup.write_potcar_file()
        kpt_setup.write_kpoints_file()
        inc_setup.write_predefined_incar_file(simulation_type=DFT_run_type, xc_functional=type_of_exchange_correlation_functional,
                                       number_of_nodes=number_of_nodes, disable_symmetry=disable_symmetry, poscar="POSCAR", material_type=electronic_structure_type,
                                       write_chgcar=write_chgcar_file, write_wavecar=write_wavecar_file)
        sub_setup.write_submit_file()

        if submit_new_jobs == bool(True):
            JobSubmission().submit_job(directory=new_dir, submit_file='submit.sh', queue_command='sbatch')

if __name__=="__main__":
    main()