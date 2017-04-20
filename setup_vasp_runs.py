__author__ = 'Ryan'

from VASP_Setup import PotcarFileSetup, KpointsFileSetup, IncarFileSetup, SubmitFileSetup
from VASP_Analyzer import DirectoryUtilities
import os
import logging
import warnings

def main():
    directory_list = DirectoryUtilities().get_downmost_directory_list()
    logging.basicConfig(filename='setup_vasp_runs.log', level=logging.INFO)

    for directory in directory_list:
        os.chdir(directory)

        # generate potcars
        if os.path.exists(directory+"/"+"POSCAR"):
            pfs = PotcarFileSetup(psp_type="PBE", poscar="POSCAR")
            kptset = KpointsFileSetup(mesh_type="Monkhorst-Pack", mesh_value="2 2 2")
            incarset = IncarFileSetup()
            submitset = SubmitFileSetup(cluster="ACI", queue="morgan", number_of_nodes=2)
            kptset.write_kpoints_file()
            incarset.write_predefined_incar_file(simulation_type="full_relax", xc_functional="GGA", number_of_nodes=2,
                                                 disable_symmetry=True)
            submitset.write_submit_file()
            pfs.write_potcar_file()
            logging.info('Successfully created input files in directory '+directory)
        else:
            warnings.warn('Could not create POTCAR file in directory '+directory)
            logging.info('Could not create POTCAR file in directory '+directory)

if __name__=='__main__':
    main()