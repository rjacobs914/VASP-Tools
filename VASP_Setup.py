__author__ = 'Ryan Jacobs'
__version__ = '2.0'
__date__ = 'Last updated March 27, 2017'
"""
VASP_Setup is a module designed to write the necessary VASP input files to conduct a DFT calculation, namely the
POSCAR, POTCAR, KPOINTS and INCAR and submission script files. There are also additional routines to aid in construction
of POSCAR files corresponding to different physical situations of interest beyond perfect unit cells, including:
supercells, surfaces, strain, and defects.

Usage note: it is necessary to set the location of the VASP pseudopotentials as the environment variable VASP_PSP_DIR.
At the Unix terminal, open your ~/.bashrc file, and set the variable by writing 'export VASP_PSP_DIR=/path_to_vasp_psp/'
"""

import math
import os
import gzip
import shutil
import subprocess
from pymatgen.matproj.rest import MPRester
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.core import Structure
from pymatgen.transformations.standard_transformations import SupercellTransformation, SubstitutionTransformation, \
    OrderDisorderedStructureTransformation,  OxidationStateDecorationTransformation
from pymatgen.transformations.defect_transformations import VacancyTransformation, InterstitialTransformation
from mpinterfaces.interface import Interface
from VASP_Analyzer import PoscarAnalyzer

class KpointsFileSetup(object):
    """Class used to create VASP KPOINTS file.
    args:
        mesh_type : (str) Type of mesh to use for kpoints.
            mesh_type = "Monkhorst-Pack" : MP-type kpt grid
            mesh_type = "Gamma" : Gamma-type kpt grid
            mesh_type = "Auto" : Automatically generate grid based on reciprocal lattice vector lengths
        mesh_values : (str) Density of kpt mesh used
            mesh_values = "4 4 4" : This is an example mesh if you're using Monkhorst-Pack or Gamma mesh.
            mesh_values = "50" : This is the kpt density if you're using Auto. Approx values are 50 for insulators,
                100 for d-band metals. Cut in half if doing HSE runs
    instance methods:
        write_kpoints_file : (None) Writes the KPOINTS file to current directory
    """

    def __init__(self, mesh_type, mesh_value, title="New kpoints file"):
        self.mesh_type = mesh_type
        self.mesh_value = mesh_value
        self.title = title

    def write_kpoints_file(self):
        if self.mesh_type == "Monkhorst-Pack" or self.mesh_type == "Gamma":
            kpoints_file = open("KPOINTS", "w")
            kpoints_file.write(str(self.title)+"\n")
            kpoints_file.write(str("0")+"\n")
            kpoints_file.write(str(self.mesh_type)+"\n")
            kpoints_file.write(str(self.mesh_value)+"\n")
            kpoints_file.write("0 0 0"+"\n")
            kpoints_file.close()
        elif self.mesh_type == "Auto":
            kpoints_file = open("KPOINTS", "w")
            kpoints_file.write(str(self.title)+"\n")
            kpoints_file.write(str("0")+"\n")
            kpoints_file.write(str(self.mesh_type)+"\n")
            kpoints_file.write(str(self.mesh_value)+"\n")
            kpoints_file.close()
        else:
            raise ValueError('mesh_type must equal either "Auto", "Gamma", or "Monkhorst-Pack"')
        return None

class PotcarFileSetup(object):
    """Class used to create a VASP POTCAR file, using elements designated in an
    existing POSCAR file.
    args:
        psp_type : (str) Type of pseudopotentials to be used in the VASP calculation
            psp_type = "PBE": pbe-type pseudopotentials. These are required if conducting phase stability analysis
                of if conducting HSE-type calculations.
            psp_type = "PW91": pw91-type pseudopotentials
        poscar : (str) The name of a poscar file. Defaults to "POSCAR"
    instance methods:
        write_potcar_file : (None) Writes a POTCAR file to current directory
    """

    def __init__(self, psp_type, poscar="POSCAR"):
        self.psp_type = psp_type
        self.poscar = poscar

    def write_potcar_file(self):
        vasp_psp_dir = os.environ["VASP_PSP_DIR"]

        if self.psp_type == "PBE" or self.psp_type == "pbe":
            psp_path = vasp_psp_dir+"/POT_GGA_PAW_PBE/"
        elif self.psp_type == "PW91" or self.psp_type == "pw91":
            psp_path = vasp_psp_dir+"/POT_GGA_PAW_PW91/"
        else:
            raise ValueError('psp_path must equal either "PBE" or "PW91"')

        if os.path.exists(os.getcwd()+"/"+"POSCAR"):
            element_names = PoscarAnalyzer(self.poscar).get_element_names()
        else:
            raise IOError('No POSCAR file exists in the working directory')

        potcar_path_list = []
        for entry in element_names:
            #Need to use element names to get list of Potcar labels (e.g., _sv, _pv, etc for each element)
            potcar_labels = self._get_potcar_labels_dict()
            potcar_type = potcar_labels[entry]
            potcar_str = "POTCAR."+potcar_type+".gz"
            potcar_path_list.append(psp_path+potcar_str)

        potcar_file = open("POTCAR", "w")
        for entry in potcar_path_list:
            potcar_temp = gzip.open(entry, "r")
            for line in potcar_temp:
                potcar_file.write(line)
        potcar_file.close()

        return None

    def _get_potcar_labels_dict(self):
        #Using POTCARs consistent with Materials Project
        if self.psp_type == "PBE" or self.psp_type == "pbe":
            potcar_labels_dict = {"Li": "Li_sv", "Na": "Na_pv", "K": "K_sv", "Cs": "Cs_sv", "Rb": "Rb_sv",
            "Be": "Be_sv", "Mg": "Mg_pv", "Ca": "Ca_sv", "Sr": "Sr_sv", "Ba": "Ba_sv", "Sc": "Sc_sv", "Y": "Y_sv",
            "Ti": "Ti_sv", "Zr": "Zr_sv", "Hf": "Hf_pv", "V": "V_sv", "Nb": "Nb_pv", "Ta": "Ta_pv", "Cr": "Cr_pv",
            "Mo": "Mo_pv", "W": "W_pv", "Mn": "Mn_pv", "Tc": "Tc_pv", "Re": "Re_pv", "Fe": "Fe_pv", "Co": "Co",
            "Ni": "Ni_pv", "Cu": "Cu_pv", "Zn": "Zn", "Ru": "Ru_pv", "Rh": "Rh_pv", "Pd": "Pd", "Ag": "Ag", "Cd": "Cd",
            "Hg": "Hg", "Au": "Au", "Ir": "Ir", "Pt": "Pt", "Os": "Os_pv", "Ga": "Ga_d", "Ge": "Ge_d", "Al": "Al",
            "As": "As", "Se": "Se", "In": "In_d", "Sn": "Sn_d", "Tl": "Tl_d", "Pb": "Pb_d", "Bi": "Bi_d",
            "Po": "Po", "At": "At_d", "La": "La", "Ce": "Ce", "Pr": "Pr_3", "Nd": "Nd_3", "Pm": "Pm_3", "Sm": "Sm_3",
            "Eu": "Eu", "Gd": "Gd", "Tb": "Tb_3", "Dy": "Dy_3", "Ho": "Ho_3", "Er": "Er_3", "Tm": "Tm_3", "Yb": "Yb",
            "Lu": "Lu_3", "B": "B", "C": "C", "N": "N", "F": "F", "Cl": "Cl", "I": "I", "Br": "Br", "O": "O", "H": "H",
            "P": "P", "Sb": "Sb", "Si": "Si", "S": "S", "Te": "Te"}
        elif self.psp_type == "PW91" or self.psp_type == "pw91":
            potcar_labels_dict = {"Li": "Li_sv", "Na": "Na_pv", "K": "K_sv", "Cs": "Cs_sv", "Rb": "Rb_sv",
            "Be": "Be_sv", "Mg": "Mg_pv", "Ca": "Ca_sv", "Sr": "Sr_sv", "Ba": "Ba_sv", "Sc": "Sc_sv", "Y": "Y_sv",
            "Ti": "Ti_sv", "Zr": "Zr_sv", "Hf": "Hf_pv", "V": "V_sv", "Nb": "Nb_pv", "Ta": "Ta_pv", "Cr": "Cr_pv",
            "Mo": "Mo_pv", "W": "W_pv", "Mn": "Mn_pv", "Tc": "Tc_pv", "Re": "Re_pv", "Fe": "Fe_pv", "Co": "Co",
            "Ni": "Ni_pv", "Cu": "Cu_pv", "Zn": "Zn", "Ru": "Ru_pv", "Rh": "Rh_pv", "Pd": "Pd", "Ag": "Ag", "Cd": "Cd",
            "Hg": "Hg", "Au": "Au", "Ir": "Ir", "Pt": "Pt", "Os": "Os_pv", "Ga": "Ga_d", "Ge": "Ge_d", "Al": "Al",
            "As": "As", "Se": "Se", "In": "In_d", "Sn": "Sn_d", "Tl": "Tl_d", "Pb": "Pb_d", "Bi": "Bi_d",
            "Po": "Po", "At": "At_d", "La": "La", "Ce": "Ce", "Pr": "Pr_3", "Nd": "Nd_3", "Pm": "Pm_3", "Sm": "Sm_3",
            "Eu": "Eu", "Gd": "Gd_3", "Tb": "Tb_2", "Dy": "Dy_3", "Ho": "Ho_3", "Er": "Er_3", "Tm": "Tm_3", "Yb": "Yb",
            "Lu": "Lu_3", "B": "B", "C": "C", "N": "N", "F": "F", "Cl": "Cl", "I": "I", "Br": "Br", "O": "O", "H": "H",
            "P": "P", "Sb": "Sb", "Si": "Si", "S": "S", "Te": "Te"}
        else:
            raise ValueError('You have not specified a valid pseudopotential type. Use psp_type="PBE" or "PW91"')

        return potcar_labels_dict

class IncarFileSetup():
    """Class used to create a VASP INCAR file.
    args:
        None. See instance methods.
    instance methods:
        write_predefined_incar_file : (None) Writes an INCAR file to the current directory using predefined INCAR tags
            args:
                simulation_type : (str) Type of VASP simulation to conduct. See below for choices.
                    run_type = "full_relax": Full volume + ions relax
                    run_type = "ions_relax": Fixed volume relax, but relaxes ions
                    run_type = "static_run": Single iteration, fixed volume and ions
                    run_type = "surface_run": Surface calculation with dipole corrections, local potential data
                XC_type : (str) Type of exchange-correlation functional used. See below for choices.
                    XC_type = "GGA": Standard GGA
                    XC_type = "GGA+U": GGA with U values for transition metals consistent with Materials Project usage
                    XC_type = "HSE": Hybrid GGA with standard 25% Hartree-Fock mixing
                number_of_nodes : (int) The number of nodes to use in the VASP calculation
                disable_symmetry : (bool) Whether to remove symmetry operations from DFT calculation
                poscar : (str) The name of a poscar file
                material_type : (str) Specify whether the material is a metal or insulator. Affects some tag choices
                    material_type = "metal" : metallic material (use for s- and p- band materials)
                    material_type = "insulator" : insulating material (use for semiconductors, band and Mott insulators, and correlated metals)
                write_chgcar : (bool) Whether to write the CHGCAR file
                write_wavecar: (bool) Whether to write the WAVECAR file
        write_custom_incar_file : (None) Writes an INCAR file to current directory based on user-provided dict of INCAR tags
            args:
                incar_dict (dict) A dict of INCAR file tags (keys) and tag values (values)
    """

    def __init__(self):
        pass

    def write_predefined_incar_file(self, simulation_type, xc_functional, number_of_nodes, cores_per_node,
                                    disable_symmetry, poscar="POSCAR", material_type="insulator", write_chgcar=False,
                                    write_wavecar=False):
        if os.path.exists(os.getcwd()+"/"+poscar):
            element_names = PoscarAnalyzer(poscar).get_element_names()
        else:
            raise IOError('No POSCAR file exists in the working directory')

        incar_dict = {"ENCUT": 500, "IBRION": 2, "ISMEAR": 0, "ISPIN": 2, "LORBIT": 11, "NEDOS": 2000, "PREC": "Accurate",
                      "NELM": 250, "TIME": 0.05, "ALGO": "All"}
        if not simulation_type == 'HSE':
            incar_dict.update(NPAR=number_of_nodes)
        if simulation_type == 'HSE':
            incar_dict.update(NPAR=number_of_nodes*cores_per_node)
        if disable_symmetry == bool(True):
            incar_dict.update(ISYM=0)

        if simulation_type == "full_relax":
            incar_dict.update(ISIF=3, NSW=200)
        elif simulation_type == "ions_relax":
            incar_dict.update(ISIF=2, NSW=200)
        elif simulation_type == "static_run":
            incar_dict.update(ISIF=2, NSW=0)
        elif simulation_type == "surface_run":
            incar_dict.update(ISIF=2, NSW=200, LDIPOL=True, IDIPOL=3, LVTOT=True, NELMDL=10)
        else:
            raise ValueError('You have not specified a valid setting for "simulation_type". "simulation_type" must be'
                             'one of: "full_relax", "ions_relax", "static_run", "surface_run')

        if write_chgcar == bool(False):
            incar_dict.update(LCHARG=False)
        elif write_chgcar == bool(True):
            incar_dict.update(LCHARG=True)
        if write_wavecar == bool(False):
            incar_dict.update(LWAVE=False)
        elif write_wavecar == bool(True):
            incar_dict.update(LWAVE=True)

        if material_type == "insulator":
            incar_dict.update(SIGMA=0.05)
        elif material_type == "metal":
            incar_dict.update(SIGMA=0.2)
        else:
            raise ValueError('"material_type" must be set to either "metal" or "insulator"')

        if xc_functional == "GGA-PW91":
            incar_dict.update(VOSKOWN=1)
        elif xc_functional == "GGA":
            pass
        elif xc_functional == "HSE":
            incar_dict.update(LHFCALC=True)
            incar_dict.update(HFSCREEN=0.2)
            incar_dict.update(PRECFOCK="Low")
            incar_dict.update(AEXX=0.25)
        elif xc_functional == "GGA+U":
            incar_dict.update(LDAU=True)
            incar_dict.update(LDAUTYPE=2)
            TM_GGAU_dict = {"V": "3.25", "Cr": "3.7", "Mn": "3.9", "Fe": "5.3", "Co": "3.32", "Ni": "6.2", "Mo": "4.38"}
            ldauj_list = []; ldaul_list = []; ldauu_list = [];
            for entry in element_names:
                if entry in TM_GGAU_dict:
                    ldauj_list.append("0")
                    ldaul_list.append("2")
                    ldauu_list.append(str(TM_GGAU_dict[entry]))
                else:
                    ldauj_list.append("0")
                    ldaul_list.append("-1")
                    ldauu_list.append("0")
            ldauj_str = ''; ldauu_str = ''; ldaul_str = ''
            for index in range(len(ldauj_list)):
                ldauj_str += str(ldauj_list[index])+" "
            for index in range(len(ldaul_list)):
                ldaul_str += str(ldaul_list[index])+" "
            for index in range(len(ldauu_list)):
                ldauu_str += str(ldauu_list[index])+" "
            incar_dict.update(LDAUJ=ldauj_str)
            incar_dict.update(LDAUU=ldauu_str)
            incar_dict.update(LDAUL=ldaul_str)
        else:
            raise ValueError('"xc_functional" must be set to one of "GGA", "GGA-PW91", "HSE", or "GGA+U"')

        incar_file = open("INCAR", "w")
        for key, value in incar_dict.items():
            incar_file.write(str(key)+" = "+str(value)+"\n")
        incar_file.close()
        return incar_dict

    def write_custom_incar_file(self, incar_dict):
        incar_file = open("INCAR", "w")
        for key, value in incar_dict.items():
            incar_file.write(str(key)+" = "+str(value)+"\n")
        incar_file.close()
        return incar_dict

class SubmitFileSetup(object):
    """Class used to create a cluster submission script for VASP runs.
    args:
        cluster : (str) The name of the computing cluster to be used for the calculation. Select from:
            cluster = "Turnbull" : run on Turnbull
            cluster = "ACI" : run on UW-Madison ACI
            cluster = "Stampede" : run on TACC-Stampede
            cluster = "Cori" : run on NERSC-Cori
            cluster = "MRSEC" : run on MRSEC cluster
        queue : (str) The name of the submission queue on the cluster of choice
            queue = "morgan.q" if cluster = "Turnbull"
            queue = "univ", "univ2", "morgan", or "morgan2" if cluster = "ACI"
            queue = "normal" if cluster = "Stampede"
            queue = "regular" if cluster = "Cori"
            queue = "all.q" if cluster = "MRSEC"
        number_of_nodes : (int) The number of nodes to be requested for the calculation
        gamma_point : (bool) Whether to use the gamma-point compiled version of VASP. This will only work if your KPOINTS
            file also uses the gamma-type mesh.
        queue_partition: (str) name of queue partition, currently for Cori only
            queue_partition = "haswell" to use haswell partition
            queue_partition = "knl" to use KNL partition
    instance methods:
        write_submit_file (None) Writes a submit.sh file to the current directory
    """
    def __init__(self, cluster, queue, number_of_nodes, gamma_point=False, queue_partition=None):
        self.cluster = cluster
        self.queue = queue
        self.number_of_nodes = number_of_nodes
        self.gamma_point = gamma_point
        self.queue_partition = queue_partition

    def write_submit_file(self):
        submit_file = open("submit.sh", "w")
        if self.cluster == "Cori" or self.cluster == "cori":
            cores_per_node = 24
            number_of_cores = int(self.number_of_nodes)*cores_per_node #24 cores per node on Cori at NERSC
            submit_file.write("#!/bin/bash"+"\n")
            submit_file.write("#SBATCH -J my_job"+"\n")
            submit_file.write("#SBATCH -p "+str(self.queue)+"\n")
            submit_file.write("#SBATCH -C "+str(self.queue_partition)+"\n")
            submit_file.write("#SBATCH -N "+str(self.number_of_nodes)+"\n")
            submit_file.write("#SBATCH -t 48:00:00"+"\n")
            submit_file.write("#SBATCH --ntasks-per-node="+str(cores_per_node)+"\n")
            submit_file.write("module load vasp/5.4.1"+"\n")
            if self.gamma_point == False:
                submit_file.write("srun -n "+str(number_of_cores)+" vasp_std > vasp.out"+"\n")
            if self.gamma_point == True:
                submit_file.write("srun -n "+str(number_of_cores)+" vasp_gam > vasp.out"+"\n")
        elif self.cluster == "stampede" or self.cluster == "Stampede":
            cores_per_node = 16
            number_of_cores = int(self.number_of_nodes)*cores_per_node #16 cores per node on Stampede
            submit_file.write("#!/bin/bash"+"\n")
            submit_file.write("#SBATCH -J my_job"+"\n")
            #submit_final.write("#SBATCH -O H.log"+"\n")
            submit_file.write("#SBATCH -p "+str(self.queue)+"\n")
            submit_file.write("#SBATCH -N "+str(self.number_of_nodes)+"\n")
            submit_file.write("#SBATCH -n "+str(number_of_cores)+"\n")
            submit_file.write("#SBATCH -t 48:00:00"+"\n")
            submit_file.write("#SBATCH --ntasks-per-node="+str(cores_per_node)+"\n")
            submit_file.write("module load python"+"\n")
            submit_file.write("module load vasp/5.4.1"+"\n")
            if self.gamma_point == False:
                submit_file.write("ibrun $TACC_VASP_DIR/bin/vasp_std > vasp.out"+"\n")
            if self.gamma_point == True:
                submit_file.write("ibrun $TACC_VASP_DIR/bin/vasp_gam > vasp.out"+"\n")
        elif self.cluster == "ACI" or self.cluster == "aci":
            cores_per_node = 16
            number_of_cores = int(self.number_of_nodes)*cores_per_node #16 cores per node on ACI
            submit_file.write("#!/bin/bash"+"\n")
            submit_file.write("#SBATCH -J my_job"+"\n")
            submit_file.write("#SBATCH -p "+str(self.queue)+"\n")
            submit_file.write("#SBATCH -N "+str(self.number_of_nodes)+"\n")
            submit_file.write("#SBATCH -n "+str(number_of_cores)+"\n")
            submit_file.write("#SBATCH -t 168:00:00"+"\n")
            submit_file.write("#SBATCH --ntasks-per-node="+str(cores_per_node)+"\n")
            submit_file.write("module load mpi/intel/mvapich2-1.9"+"\n")
            submit_file.write("module load compile/intel"+"\n")
            if self.gamma_point == False:
                submit_file.write("/usr/mpi/intel/mvapich2-1.9/bin/mpiexec -n "+str(number_of_cores)+" /home/groups/CMG/bin/vasp.5.4.1/vasp_std > vasp.out"+"\n")
            if self.gamma_point == True:
                submit_file.write("/usr/mpi/intel/mvapich2-1.9/bin/mpiexec -n "+str(number_of_cores)+" /home/groups/CMG/bin/vasp.5.4.1/vasp_gam > vasp.out"+"\n")
        elif self.cluster == "Turnbull" or self.cluster == "turnbull":
            cores_per_node = 12
            number_of_cores = int(self.number_of_nodes)*cores_per_node #12 cores per node on Turnbull
            submit_file.write("#!/bin/bash"+"\n")
            submit_file.write("#$ -N my_job"+"\n")
            submit_file.write("#$ -q "+str(self.queue)+"\n")
            submit_file.write("#$ -pe mvapich2 "+str(number_of_cores)+"\n")
            submit_file.write("#$ -l h_rt=168:00:00"+"\n")
            submit_file.write("#$ -cwd"+"\n")
            submit_file.write("#$ -o $JOB_NAME.o$JOB_ID"+"\n")
            submit_file.write("#$ -e $JOB_NAME.e$JOB_ID"+"\n")
            submit_file.write("source /opt/intel/composer_xe_2011_sp1.7.256/mkl/bin/mklvars.sh intel64"+"\n")
            if self.gamma_point == False:
                submit_file.write("/usr/local/mvapich2/intel/1.8.1/bin/mpiexec -n "+str(number_of_cores)+" /share/apps/bin/vasp/5.4.1/vasp_std > vasp.out"+"\n")
            if self.gamma_point == True:
                submit_file.write("/usr/local/mvapich2/intel/1.8.1/bin/mpiexec -n "+str(number_of_cores)+" /share/apps/bin/vasp/5.4.1/vasp_gma > vasp.out"+"\n")

        elif self.cluster == "MRSEC" or self.cluster == "mrsec":
            cores_per_node = 12
            number_of_cores = int(self.number_of_nodes)*cores_per_node #12 cores per node on MRSEC
            submit_file.write("#!/bin/bash"+"\n")
            submit_file.write("#$ -S /bin/bash"+"\n")
            submit_file.write("#$ -N my_job"+"\n")
            submit_file.write("#$ -cwd"+"\n")
            submit_file.write("#$ -V"+"\n")
            submit_file.write("#$ -q "+str(self.queue)+"\n")
            submit_file.write("#$ -pe mpich "+str(number_of_cores)+"\n")
            submit_file.write("#$ -l h_rt=168:00:00"+"\n")
            submit_file.write("#$ -o $JOB_NAME.o$JOB_ID"+"\n")
            submit_file.write("#$ -e $JOB_NAME.e$JOB_ID"+"\n")
            submit_file.write("cat $TMPDIR/machines | uniq>mycluster"+"\n")
            submit_file.write("MPI_HOME=/usr/local/mvapich2/intel/1.8.1/bin"+"\n")
            if self.gamma_point == False:
                submit_file.write("$MPI_HOME/mpiexec -n $NSLOTS /share/apps/vasp/vasp-5.3.3/std/vasp > vasp.out"+"\n")
            if self.gamma_point == True:
                submit_file.write("$MPI_HOME/mpiexec -n $NSLOTS /share/apps/vasp/vasp-5.3.3/gma/vasp > vasp.out"+"\n")
        else:
            raise ValueError('"cluster" was not defined properly. Pick one of "MRSEC", "Turnbull", "ACI", "Cori", or "Stampede')

        submit_file.close()
        return None

class PoscarFileSetup(object):
    """Class used to create a VASP POSCAR file by fetching structure data from Materials Project
    args:
        material_composition : (str) The chemical formula of the material to make a POSCAR file
        mapi_key : (str) Your Materials API key from MaterialsProject.org
        get_only_most_stable_structure : (bool) Whether to only return most stable structure when querying MaterialsProject
    instance methods:
        write_poscar_file : (list of poscar objects) Writes POSCAR file(s) to current directory
    """
    def __init__(self, material_composition, mapi_key, get_only_most_stable_structure=True):
        self.material_composition = material_composition
        self.mapi_key = mapi_key
        self.get_only_most_stable_structure = get_only_most_stable_structure

    def write_poscar_file(self, structure_list=None):
        if structure_list is None:
            structure_list, structure_data_list = self._get_structure_from_materials_project()

        stable_materialid = self._get_most_stable_materialid(structure_data_list=structure_data_list)

        poscar_list = []
        if self.get_only_most_stable_structure == bool(False):
            for index, structure in enumerate(structure_list):
                poscar = Poscar(structure=structure)
                poscar_file = open("POSCAR"+"_"+str(structure_data_list[index]["material_id"]), "w")
                poscar_file.write(str(poscar))
                poscar_file.close()
                poscar_list.append(poscar)
        elif self.get_only_most_stable_structure == bool(True):
            for index, structure in enumerate(structure_list):
                if structure_data_list[index]["material_id"] == stable_materialid:
                    poscar = Poscar(structure=structure)
                    poscar_file = open("POSCAR", "w")
                    poscar_file.write(str(poscar))
                    poscar_file.close()
                    poscar_list.append(poscar)
        else:
            raise ValueError('You must set "get_only_most_stable_structure" to either True or False (as bool)')

        return poscar_list

    def _get_most_stable_materialid(self, structure_data_list):
        lowest_energy = 1000
        for structure_dict in structure_data_list:
            if structure_dict["e_above_hull"] < lowest_energy:
                stable_materialid = structure_dict["material_id"]
                lowest_energy = structure_dict["e_above_hull"]
        return stable_materialid

    def _get_structure_from_materials_project(self):
        material_ids = []
        mprester = MPRester(self.mapi_key)
        structure_data_list = self._get_data_from_materials_project()
        structure_list = []; structure_conventional_list = []
        for structure_dict in structure_data_list:
            structure = mprester.get_structure_by_material_id(material_id=str(structure_dict["material_id"]))
            structure_list.append(structure)
            material_ids.append(str(structure_dict["material_id"]))

        # Convert the primitive cell to a conventional cell for easier handling
        for structure in structure_list:
            sa = SpacegroupAnalyzer(structure=structure)
            structure_conventional = sa.get_conventional_standard_structure()
            structure_conventional_list.append(structure_conventional)
        return structure_conventional_list, structure_data_list

    def _get_data_from_materials_project(self):
        mprester = MPRester(self.mapi_key)
        structure_data_list = []
        structure_data = mprester.get_data(chemsys_formula_id=self.material_composition)
        for structure_data_dict in structure_data:
            structure_data_dict_condensed = {}
            structure_data_dict_condensed["e_above_hull"]=structure_data_dict["e_above_hull"]
            structure_data_dict_condensed["spacegroup_symbol"]=structure_data_dict["spacegroup"]["symbol"]
            structure_data_dict_condensed["spacegroup_number"]=structure_data_dict["spacegroup"]["number"]
            structure_data_dict_condensed["pretty_formula"]=structure_data_dict["pretty_formula"]
            structure_data_dict_condensed["material_id"]=structure_data_dict["material_id"]
            structure_data_list.append(structure_data_dict_condensed)
        return structure_data_list

class PoscarFileModifier(object):
    """Class used to perform various useful structural modifications to a basic POSCAR file
    args:
        poscar : (str) Name of a poscar file in current directory
    instance methods:
        make_poscar_supercell : (poscar object) Creates a supercell from an existing POSCAR file
            args:
                scaling_matrix : (list) A list of three scaling vectors for supercell dimensions, e.g. [[2, 0, 0],[0, 2, 0],[0, 0, 3]]
        make_poscar_rescaled : (poscar object) Function that rescales the supercell dimensions to incorporate strain, and to also simply rescale the lattice
            while keeping all atoms fixed at their relative direct coordinates. Useful for creating strained materials,
            rescaling vacuum dimensions for surfaces or nanoparticles.
            args:
                cell_rescale_dict : (dict) A dictionary of cell directions and cell lengths (in Angstroms) pairs to rescale
                    the cell by, e.g. {"x": 6, "z": 40}. The cell directions must be one of "x", "y" and "z" and the cell
                    lengths must be in units of Angstroms.
                is_strained : (bool) Whether to treat the rescaling as a strain (i.e., atom positions are scaled with the
                    new cell lengths) (is_strained = True) or whether to simply rescale the cell length and keep atom
                    positions fixed (is_strained = False). Note that this type of rescaling is best for changing the size of
                    surface cells along the vacuum direction, and some unphysical structure may result when rescaling in
                    the other directions.
        make_poscar_with_dopants : (poscar object) This function is used to create an alloyed POSCAR, and more generally
            to easily change the composition of a POSCAR
            args:
                doping_type : (str) Designates rather to do single or multi-site alloying
                    doping_type = "Full", alloying where one species is converted completely to another species. Can be done on multiple sites
                    doping_type = "Partial", alloying where one site is converted to multiple species with partial occupancies. Can be done on multiple sites
                doping_scheme : (dict) A dictionary of species to alloy, see examples below:
                    if doping_type = "Full", doping_scheme={"Li": "Na", "Cl": "F"} will replace all Li with Na and all Cl with F.
                    if doping_type = "Partial", doping_scheme={"Li": {"Na": 0.75, "K": 0.25}} will replace Li with 75% Na and 25% K. Note that for determining the ordered
                        structure with different occupancies, this routine can be computationally expensive.
                oxidation states : (dict) A dictionary of oxidation states,
                    if doping_type = "Partial", {"Li":1, "Na":1, "K":1 "O":-2}. Need to supply oxidation states of all elements involved
        make_poscar_with_vacancy : (list of poscar objects) This function is used to create POSCAR files with each symmetry distinct vacancy, up to a specified supercell size.
            args:
                vacancy_species: (list) A list of elemental species to create vacancies, e.g. ["Li", "Co"]
                max_supercell_size : (list of list) A matrix of the maximum supercell size to generate defected structures.
                    e.g., [[2, 0, 0], [0, 2, 0], [0, 0, 2]]
        make_poscar_slab : (poscar object) This function is used to create a surface slab from an existing POSCAR file
            args:
                surface_hkl: (list) The miller indices with which to create the surface slab, e.g. [1, 1, 0]
                minimum_thickness: (int) The minimum thickness of the resulting supercell slab, in the direction of the specified
                    miller indices hkl.
                minimum_vacuum: (int) The minimum thickness of the vacuum region, in the direction of the specified miller indices hkl.
                shift_slab_to_bottom: (bool) Whether to shift the atom coordinates of the slab such that the slab is
                    positioned at the bottom of the supercell. This can sometimes aid in surface convergence with dipole
                    corrections.
        trim_poscar_slab : (poscar object) This function is used to trim the top or bottom layer of atoms from a surface slab.
            Useful for isolating the desired surface termination if multiple terminations of a single hkl direction are possible.
            Note to users: this routine may fail for materials containing surfaces where the atom positions vary in the
            z-direction by more than the value of atom_position_tolerance. You may need to experiment with the value of this
            parameter to get the desired results. Use carefully!
            args:
                surface_trimming_dict: (dict) A dictionary written with keys of "Top" and "Bottom" and values
                    equal to number of surface layers to trim on each surface. For example, {"Top": 2, "Bottom": 1} will remove the
                    top two surface layers and bottom one surface layers from the POSCAR.
                atom_position_tolerance: (float) cutoff distance
        add_selective_dynamics : (poscar object) This function is used to add selective dynamics tags of 'T T T' and
            'F F F' to surface atoms. Generally, one wants to relax the top and bottom couple layers, and freeze the
            remaining atoms in the slab to be bulk-like.
            args:
                number_of_layers : (int) The number of layers from the surface to designate as 'T T T', i.e. the number
                    of layers you want to relax.
                top_and_bottom : (bool) Whether or not to relax both the top and bottom surfaces. Defaults to True. If set
                    to False, only the top layer is relaxed.
                layer_thickness : (float) Surface layer thickness, in angstroms. A typical value is 2, which is set as default.
    """
    def __init__(self, poscar="POSCAR"):
        self.poscar = poscar
        shutil.copy(os.getcwd()+"/"+self.poscar, os.getcwd()+"/"+self.poscar+"_original")

    def make_poscar_supercell(self, scaling_matrix):
        poscar_placeholder = "POSCAR_placeholder"
        structure = Structure.from_file(self.poscar)
        structure_scaled = SupercellTransformation(scaling_matrix).apply_transformation(structure)
        poscar = Poscar(structure_scaled)
        poscar.write_file(poscar_placeholder)
        shutil.move(os.getcwd()+"/"+poscar_placeholder, os.getcwd()+"/"+self.poscar)
        return poscar

    def make_poscar_rescaled(self, cell_rescale_dict, is_strained=False):
        psr = PoscarAnalyzer(poscar=self.poscar)
        lattice_parameters = psr.get_lattice_parameters
        atom_positions = psr.get_atom_positions_direct
        element_names = psr.get_element_names
        atom_amounts = psr.get_atom_amounts

        for key in cell_rescale_dict.keys():
            if key == "x":
                new_cell_length = float(cell_rescale_dict[key])
                old_cell_length = math.sqrt((lattice_parameters[0][0])**2+(lattice_parameters[0][1])**2+(lattice_parameters[0][2])**2)
                lattice_parameters[0][0] = lattice_parameters[0][0]*(new_cell_length/old_cell_length)
                if is_strained == bool(False):
                    for index in range(len(atom_positions[-1])):
                        atom_positions[0][index] *= atom_positions[0][index]*(old_cell_length/new_cell_length)

            if key == "y":
                new_cell_length = float(cell_rescale_dict[key])
                old_cell_length = math.sqrt((lattice_parameters[1][0])**2+(lattice_parameters[1][1])**2+(lattice_parameters[1][2])**2)
                lattice_parameters[1][1] = lattice_parameters[1][1]*(new_cell_length/old_cell_length)
                if is_strained == bool(False):
                    for index in range(len(atom_positions[-1])):
                        atom_positions[1][index] = atom_positions[1][index]*(old_cell_length/new_cell_length)

            if key == "z":
                new_cell_length = float(cell_rescale_dict[key])
                old_cell_length = math.sqrt((lattice_parameters[2][0])**2+(lattice_parameters[2][1])**2+(lattice_parameters[2][2])**2)
                lattice_parameters[2][2] = lattice_parameters[2][2]*(new_cell_length/old_cell_length)
                if is_strained == bool(False):
                    for index in range(len(atom_positions[-1])):
                        atom_positions[2][index] = atom_positions[2][index]*(old_cell_length/new_cell_length)

        poscar_original = open(self.poscar, "r")
        poscar_placeholder = "POSCAR_placeholder"
        poscar_scaled = open(poscar_placeholder, "w")
        line_count = 0
        for line in poscar_original:
            if line_count < 2:
                poscar_scaled.write(line)
                line_count += 1
        if line_count == 2:
            poscar_scaled.write(str(lattice_parameters[0][0])+" "+str(lattice_parameters[0][1])+" "+str(lattice_parameters[0][2])+"\n")
            line_count += 1
        if line_count == 3:
            poscar_scaled.write(str(lattice_parameters[1][0])+" "+str(lattice_parameters[1][1])+" "+str(lattice_parameters[1][2])+"\n")
            line_count += 1
        if line_count == 4:
            poscar_scaled.write(str(lattice_parameters[2][0])+" "+str(lattice_parameters[2][1])+" "+str(lattice_parameters[2][2])+"\n")
            line_count += 1
        if line_count == 5:
            poscar_scaled.write(str(element_names[0])+" "+str(element_names[1])+" "+str(element_names[2])+"\n")
            line_count += 1
        if line_count == 6:
            poscar_scaled.write(str(atom_amounts[0])+" "+str(atom_amounts[1])+" "+str(atom_amounts[2])+"\n")
            line_count += 1
        if line_count == 7:
            poscar_scaled.write("Selective Dynamics"+"\n")
            line_count += 1
        if line_count == 8:
            poscar_scaled.write("Direct"+"\n")
            line_count += 1
        if line_count > 8:
            for index in range(len(atom_positions[-1])):
                poscar_scaled.write(str(atom_positions[0][index])+" "+str(atom_positions[1][index])+" "+str(atom_positions[2][index])+"\n")
                line_count += 1

        poscar_original.close()
        poscar_scaled.close()
        shutil.move(os.getcwd() + "/" + poscar_placeholder, os.getcwd() + "/" + self.poscar)
        structure = Structure.from_file(self.poscar)
        poscar = Poscar(structure=structure)
        return poscar

    def make_poscar_with_dopants(self, doping_type, doping_scheme, oxidation_states=None):
        structure = Structure.from_file(self.poscar)
        poscar_placeholder = "POSCAR_placeholder"
        if doping_type == "Full":
            structure_alloyed = SubstitutionTransformation(species_map=doping_scheme).apply_transformation(structure=structure)
            poscar = Poscar(structure=structure_alloyed)
            poscar.write_file(poscar_placeholder)
        elif doping_type == "Partial":
            structure_alloyed = SubstitutionTransformation(species_map=doping_scheme).apply_transformation(structure=structure)
            structure_alloyed_decorated = OxidationStateDecorationTransformation(oxidation_states=oxidation_states).apply_transformation(structure=structure_alloyed)
            structure_alloyed_decorated_ordered = OrderDisorderedStructureTransformation().apply_transformation(structure=structure_alloyed_decorated)
            poscar = Poscar(structure=structure_alloyed_decorated_ordered)
            poscar.write_file(poscar_placeholder)
        else:
            raise ValueError('"doping_type" must be set to either "Full" or "Partial"')

        shutil.move(os.getcwd() + "/" + poscar_placeholder, os.getcwd() + "/" + self.poscar)
        return poscar

    def make_poscar_with_vacancy(self, vacancy_species, max_supercell_size):
        element_names = PoscarAnalyzer(poscar=self.poscar).get_element_names
        count = 1
        structure = Structure.from_file(self.poscar)
        poscar_list = []

        # Check that vacancy species is an element in the POSCAR
        for element in vacancy_species:
            if element not in element_names:
                raise ValueError('An entry in "vacancy_species" is not present in the POSCAR, cannot create vacancy structures')

        for element in vacancy_species:
            vac = VacancyTransformation(supercell_dim=max_supercell_size, species=element)
            structure_with_vacancies = vac.apply_transformation(structure=structure, return_ranked_list=True)
            poscar = Poscar(structure=structure_with_vacancies[0]["structure"])
            poscar_list.append(poscar)
            poscar.write_file("POSCAR"+"_"+str(element)+str(count))
            count += 1
        return poscar_list

    def make_poscar_slab(self, surface_hkl, minimum_thickness, minimum_vacuum=20, shift_slab_to_bottom=True):
        bulk_structure = Structure.from_file(self.poscar)
        interface = Interface(bulk_structure, hkl=surface_hkl, min_thick=minimum_thickness, min_vac=minimum_vacuum,
                              primitive=False, from_ase=True)
        interface.create_interface()
        interface.sort()
        poscar = Poscar(interface)
        poscar_placeholder = "POSCAR_placeholder"

        if shift_slab_to_bottom == bool(False):
            poscar.write_file(poscar_placeholder)
            shutil.move(os.getcwd() + "/" + poscar_placeholder, os.getcwd() + "/" + self.poscar)
        elif shift_slab_to_bottom == bool(True):
            poscar_slab_shifted = "POSCAR_slab_shifted"
            poscar.write_file(poscar_placeholder)

            psr = PoscarAnalyzer(poscar=poscar_placeholder)
            atom_coordinates = psr.get_atom_positions_direct
            total_atoms = psr.get_total_atoms
            atom_coord_min = min(atom_coordinates[2])
            for index in range(len(atom_coordinates[2])):
                atom_coordinates[2][index] -= atom_coord_min

            poscar_placeholder_file = open(poscar_placeholder, "r")
            poscar_slab_shifted_file = open(poscar_slab_shifted, "w")
            line_count = 0
            for line in poscar_placeholder_file:
                if line_count < 8:
                    poscar_slab_shifted_file.write(line)
                    line_count += 1
                if line_count >= 8 and line_count < 8+total_atoms:
                    for index in range(int(total_atoms)):
                        poscar_slab_shifted_file.write(str(atom_coordinates[0][index])+" "+str(atom_coordinates[1][index])+" "+str(atom_coordinates[2][index])+"\n")
                        line_count += 1

            poscar_placeholder_file.close()
            poscar_slab_shifted_file.close()
            shutil.move(os.getcwd() + "/" + poscar_slab_shifted, os.getcwd() + "/" + self.poscar)
            os.remove(os.getcwd()+"/"+poscar_placeholder)
            structure = Structure.from_file(self.poscar)
            poscar = Poscar(structure=structure)

        return poscar

    def trim_poscar_slab(self, surface_trimming_dict, atom_position_tolerance = 0.01):
        # Loop through top and bottom components of dict, separately removing the top and bottom layers
        for key in surface_trimming_dict.keys():
            surface_layer_count = 0
            # Only perform surface removal if number of layers to remove is greater than zero
            if surface_trimming_dict[key] > 0:
                # Keep removing layers until the desired number of layers has been reached
                while surface_layer_count < surface_trimming_dict[key]:
                    poscar_slab = open(self.poscar, "r")
                    element_names = PoscarAnalyzer(poscar=self.poscar).get_element_names
                    list_z_coords = []; list_z_coords_unique = []; trimmed_slab_coords = [];
                    # Specify where coordinates start in POSCAR, whether SD has been specified or not
                    line_where_coords_start = 8
                    for line in poscar_slab:
                        if "Selective Dynamics" in line or "selective dynamics" in line:
                            line_where_coords_start = 9
                    poscar_slab.close()

                    line_count = 0
                    poscar_slab = open(self.poscar, "r")
                    for line in poscar_slab:
                        line_count += 1
                        if line_count > line_where_coords_start: #Needs to be 9 if Selective Dynamics present, else is equal to 8
                            list_z_coords.append(line.split()[2])
                    poscar_slab.close()

                    for line in list_z_coords:
                        if line not in list_z_coords_unique:
                            list_z_coords_unique.append(line)

                    if key == "Top":
                        layer_to_remove = float(max(list_z_coords_unique))
                    if key == "Bottom":
                        layer_to_remove = float(min(list_z_coords_unique))

                    poscar_slab = open(self.poscar, "r")
                    line_count = 0
                    for line in poscar_slab:
                        line_count += 1
                        # Generate the list of atom coordinates of atoms to KEEP
                        if line_count > line_where_coords_start: #8 if no SD, 9 if SD
                            if layer_to_remove + atom_position_tolerance < float(line.split()[2]):
                                trimmed_slab_coords.append(line)
                            if float(line.split()[2]) < layer_to_remove - atom_position_tolerance:
                                trimmed_slab_coords.append(line)
                    poscar_slab.close()

                    poscar_placeholder = "POSCAR_placeholder"
                    poscar_slab_trimmed = open(poscar_placeholder, "w")
                    atom_count_list = []
                    for entry in element_names:
                        atom_count = 0
                        for line in trimmed_slab_coords:
                            if entry in line:
                                atom_count += 1
                        atom_count_list.append(atom_count)
                    print atom_count_list

                    poscar_slab = open(self.poscar, "r")
                    line_count = 0
                    for line in poscar_slab:
                        if line_count < 6:
                            poscar_slab_trimmed.write(line)
                            line_count += 1
                    if line_count == 6:
                        for entry in atom_count_list:
                            poscar_slab_trimmed.write(str(entry)+" ")
                        poscar_slab_trimmed.write("\n")
                        line_count += 1
                    if line_count == 7:
                        poscar_slab_trimmed.write("direct"+"\n")
                        line_count += 1
                    if line_count > 7:
                        for entry in trimmed_slab_coords:
                            poscar_slab_trimmed.write(str(entry))
                            line_count += 1

                    poscar_slab_trimmed.close()
                    poscar_slab.close()
                    shutil.move(os.getcwd() + "/" + poscar_placeholder, os.getcwd() + "/" + self.poscar)
                    surface_layer_count += 1
        structure = Structure.from_file(self.poscar)
        poscar = Poscar(structure=structure)
        return poscar

    def add_selective_dynamics(self, number_of_layers, top_and_bottom = True, layer_thickness=2.0):
        #Typically a z_tolerance equal to 0.5*layer_spacing - delta will make a single layer marked "T T T".
        pa = PoscarAnalyzer(poscar=self.poscar)
        lattice_parameters = pa.get_lattice_parameters
        # Normalize layer thickness by c lattice parameter
        layer_thickness = layer_thickness/math.sqrt(lattice_parameters[2][0]**2+lattice_parameters[2][1]**2+lattice_parameters[2][2]**2)
        #delta = 0.25*layer_thickness
        #z_tolerance = (0.5*layer_thickness) + (layer_thickness*(number_of_layers-1)) - delta
        z_tolerance = (0.5*layer_thickness) + (layer_thickness*(number_of_layers-1))
        atom_positions = pa.get_atom_positions_direct

        # Create selective dynamics selection based on z-axis atom position
        z_coords = []
        for entry in atom_positions[-1]:
            if entry not in z_coords:
                z_coords.append(entry)
        zcoords_min = min(z_coords)
        zcoords_max = max(z_coords)

        count_true = 0; count_false = 0; selective_dynamics = [];
        for entry in atom_positions[-1]:
            if top_and_bottom == bool(True):
                if entry < zcoords_max + z_tolerance and entry > zcoords_max - z_tolerance:
                    selective_dynamics.append('T T T')
                    count_true += 1
                if entry > zcoords_min - z_tolerance and entry < zcoords_min + z_tolerance:
                    selective_dynamics.append('T T T')
                    count_true += 1
                if entry < zcoords_max - z_tolerance and entry > zcoords_min + z_tolerance:
                    selective_dynamics.append('F F F')
                    count_false += 1
            elif top_and_bottom == bool(False):
                if entry < zcoords_max + z_tolerance and entry > zcoords_max - z_tolerance:
                    selective_dynamics.append('T T T')
                    count_true += 1
                else:
                    selective_dynamics.append('F F F')
                    count_false += 1
            else:
                raise ValueError('"top_and_bottom" must be either True or False (as bool)')

        #print "Made %i entries True and %i entries False for selective dynamics!" % (count_true, count_false)

        poscar_nosd = open(self.poscar, "r")
        poscar_placeholder = "POSCAR_placeholder"
        poscar_sd = open(poscar_placeholder, "w")
        line_count = 0
        for line in poscar_nosd:
            if line_count < 7:
                poscar_sd.write(line)
                line_count += 1
        if line_count == 7:
            poscar_sd.write("Selective Dynamics"+"\n")
            line_count += 1
        if line_count == 8:
            poscar_sd.write("Direct"+"\n")
            line_count += 1
        if line_count > 8:
            for index in range(len(atom_positions[-1])):
                poscar_sd.write(str(atom_positions[0][index])+" "+str(atom_positions[1][index])+" "+str(atom_positions[2][index])+" "+str(selective_dynamics[index])+"\n")
                line_count += 1

        poscar_nosd.close()
        poscar_sd.close()
        shutil.move(os.getcwd() + "/" + str(poscar_sd), os.getcwd() + "/" + self.poscar)
        structure = Structure.from_file(self.poscar)
        poscar = Poscar(structure=structure)
        return poscar

class JobSubmission():
    """Class used to submit vasp jobs
    args: None
    class methods:
        submit_job : (None) Submit a VASP job from current directory
    """
    @classmethod
    def submit_job(cls, directory, submit_file="submit.sh", queue_command="sbatch"):
        os.chdir(directory)
        subprocess.Popen([queue_command, submit_file]).communicate()
        return None
