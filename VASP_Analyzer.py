__author__ = 'Ryan Jacobs'
__version__ = '2.0'
__date__ = 'Last updated March 27, 2017'
"""
VASP_Analyzer is a module designed to handle basic analysis of VASP input and output files. It also contains some
miscellaneous job utilities designed to gather completed job data and assess the convergence of jobs.

Most users will not need to call any of the functions in this module.
Instead, the VASP_PostProcessing and VASP_Setup modules draw upon the functions included here.
"""

import re
import numpy as np
import os
import xlsxwriter
import time
import getpass
import subprocess
import logging
import shutil
from pymatgen.io.vasp.inputs import Poscar

class PoscarAnalyzer(object):
    """This class is used to obtain basic, useful quantities from the POSCAR file.
    args:
        poscar: (str), the name of a POSCAR file

    Example usage:
        pa = PoscarAnalyzer(poscar="POSCAR")
        total_atoms = pa.get_total_atoms
    """
    def __init__(self, poscar="POSCAR"):
        self.poscar = poscar

    def get_element_names(self):
        poscar = open(self.poscar, "r")
        poscar_data = poscar.readlines()
        element_names = re.findall(r"[\w']+", poscar_data[5])
        poscar.close()
        return element_names

    def get_atom_amounts(self):
        atom_amounts = []
        poscar = open(self.poscar, "r")
        poscar_data = poscar.readlines()
        atoms = re.findall(r"[\w']+", poscar_data[6])
        for entry in atoms:
            atom_amounts.append(int(entry))
        poscar.close()
        return atom_amounts

    def get_total_atoms(self):
        total_atoms = 0
        atom_amounts = self.get_atom_amounts
        for index in range(len(atom_amounts)):
            total_atoms += float(atom_amounts[index])
        return total_atoms

    def get_element_composition(self):
        element_composition = []
        total_atoms = self.get_total_atoms
        atom_amounts = self.get_atom_amounts
        for index in range(len(atom_amounts)):
            element_composition.append(float(atom_amounts[index]/total_atoms))
        return element_composition

    def get_lattice_parameters(self):
        lattice_parameters = []; lattice_parameters_split = []
        poscar = open(self.poscar, "r")
        poscar_data = poscar.readlines()
        for index in range(len(poscar_data)):
            if index > 1 and index < 5:
                lattice_parameters.append(poscar_data[index])
        for index in range(len(lattice_parameters)):
            lattice_parameters_split.append(lattice_parameters[index].split())
        for index in range(len(lattice_parameters_split)):
            if index == 0:
                lattice_parameters_xcolumn = np.array([float(lattice_parameters_split[index][index]), float(lattice_parameters_split[index+1][index]), float(lattice_parameters_split[index+2][index])])
            if index == 1:
                lattice_parameters_ycolumn = np.array([float(lattice_parameters_split[index-1][index]), float(lattice_parameters_split[index][index]), float(lattice_parameters_split[index+1][index])])
            if index == 2:
                lattice_parameters_zcolumn = np.array([float(lattice_parameters_split[index-2][index]), float(lattice_parameters_split[index-1][index]), float(lattice_parameters_split[index][index])])
        lattice_parameters = np.array([lattice_parameters_xcolumn, lattice_parameters_ycolumn, lattice_parameters_zcolumn])
        poscar.close()
        return lattice_parameters

    def get_cell_volume(self):
        lattice_parameters = self.get_lattice_parameters
        lattice_parameters_x = lattice_parameters[0][0]
        lattice_parameters_y = lattice_parameters[1][1]
        lattice_parameters_z = lattice_parameters[2][2]
        cell_volume = lattice_parameters_x*lattice_parameters_y*lattice_parameters_z
        return cell_volume

    def get_atom_positions_direct(self):
        atom_positions = []; atom_positions_split = []; atom_positions_split_float_x = []; atom_positions_split_float_y = [];
        atom_positions_split_float_z = []; atom_positions_array = np.empty((0,1)); atom_positions_array_x = np.empty((0,1)); atom_positions_array_y = np.empty((0,1))
        atom_positions_array_z = np.empty((0,1));

        poscar = open(self.poscar, "r")
        poscar_data = poscar.readlines()
        total_atoms = self.get_total_atoms

        atom_position_start = 7
        atom_position_end = 8
        for line in poscar_data:
            if "Selective dynamics" in line or "Selective Dynamics" in line:
                atom_position_start = 8
                atom_position_end = 9

        for index in range(len(poscar_data)):
            if index > atom_position_start and index < total_atoms+atom_position_end:
                atom_positions.append(poscar_data[index])

        for index in range(len(atom_positions)):
            atom_positions_split.append(atom_positions[index].split())

        for index in range(len(atom_positions_split)):
            atom_positions_split_float_x.append(float(atom_positions_split[index][0]))
            atom_positions_split_float_y.append(float(atom_positions_split[index][1]))
            atom_positions_split_float_z.append(float(atom_positions_split[index][2]))

        for index in range(len(atom_positions_split)):
            atom_positions_array_x = np.insert(atom_positions_array_x, index, atom_positions_split_float_x[index])
            atom_positions_array_y = np.insert(atom_positions_array_y, index, atom_positions_split_float_y[index])
            atom_positions_array_z = np.insert(atom_positions_array_z, index, atom_positions_split_float_z[index])
            atom_positions_array = np.array([atom_positions_array_x, atom_positions_array_y, atom_positions_array_z])

        poscar.close()
        return atom_positions_array

    def get_atom_positions_cartesian(self):
        x_coords_cart = []; y_coords_cart = []; z_coords_cart = [];
        atom_positions_array_x = np.empty((0,1)); atom_positions_array_y = np.empty((0,1)); atom_positions_array_z = np.empty((0,1));

        atom_positions_array = self.get_atom_positions_direct
        lattice_parameters = self.get_lattice_parameters

        for index in range(len(atom_positions_array[0])):
            x_coords_cart.append(atom_positions_array[0][index]*lattice_parameters[0][0]+atom_positions_array[1][index]*lattice_parameters[0][1]+atom_positions_array[2][index]*lattice_parameters[0][2])
            y_coords_cart.append(atom_positions_array[0][index]*lattice_parameters[1][0]+atom_positions_array[1][index]*lattice_parameters[1][1]+atom_positions_array[2][index]*lattice_parameters[1][2])
            z_coords_cart.append(atom_positions_array[0][index]*lattice_parameters[2][0]+atom_positions_array[1][index]*lattice_parameters[2][1]+atom_positions_array[2][index]*lattice_parameters[2][2])

        for index in range(len(x_coords_cart)):
            atom_positions_array_x = np.insert(atom_positions_array_x, index, x_coords_cart[index])
            atom_positions_array_y = np.insert(atom_positions_array_y, index, y_coords_cart[index])
            atom_positions_array_z = np.insert(atom_positions_array_z, index, z_coords_cart[index])
            atom_positions_array = np.array([atom_positions_array_x, atom_positions_array_y, atom_positions_array_z])

        return atom_positions_array

    def get_atom_distances(self, write_to_file=False):
        # Make poscar into pymatgen structure object
        structure = Poscar.from_file(self.poscar)
        atom_distances = structure.distance_matrix
        if write_to_file == bool(True):
            distance_matrix_file = open("distance_matrix.txt", "w")
            for index in range(len(atom_distances)):
                distance_matrix_file.write(atom_distances[index])
            distance_matrix_file.close()

        return atom_distances

class OutcarAnalyzer(object):
    """This class is used to obtain basic, useful quantities from the OUTCAR file.
    args:
        outcar: (str), the name of an OUTCAR file

    Example usage:
        oa = OutcarAnalyzer(outcar="OUTCAR")
        fermi = oa.get_fermi_energy
    """
    def __init__(self, outcar="OUTCAR"):
        self.outcar = outcar

    def get_fermi_energy(self):
        outcar = open(str(self.outcar), "r")
        outcar_data = outcar.readlines()
        Efermi_list = []
        for entry in outcar_data:
            if 'E-fermi' in entry:
                Efermi_list.append(entry)

        Efermi_datalist = []
        for entry in Efermi_list:
            Efermi_data = re.findall(r"[-+]?\d*\.\d+|\d+",entry)
            Efermi_datalist.append(Efermi_data)

        fermi_energy = float(Efermi_datalist[-1][0])
        return fermi_energy

    def get_dipole(self):
        dipole_list = []
        outcar = open(str(self.outcar), "r")
        outcar_data = outcar.readlines()

        for entry in outcar_data:
            if 'LDIPOL =      F' in entry:
                raise ValueError("LDIPOL is set to F, no dipole value for this job")
            if 'dipolmoment' in entry:
                dipole_list.append(entry)

        dipole_datalist = []
        for entry in dipole_list:
            dipole_data = re.findall(r"[-+]?\d*\.\d+|\d+",entry)
            dipole_datalist.append(dipole_data)

        dipole = float(dipole_datalist[-1][-1])
        return dipole

class OszicarAnalyzer(object):
    """
    This class is used to obtain basic, useful quantities from the OSZICAR file.
    args:
        oszicar: (str), the name of an OSZICAR file
        poscar: (str), a POSCAR file

    Example usage:
        oza = OszicarAnalyzer(oszicar="OSZICAR", poscar="POSCAR")
        energy = oza.get_energy
    """
    def __init__(self, oszicar="OSZICAR", poscar="POSCAR"):
        self.oszicar = oszicar
        self.poscar = poscar

    def get_energy(self):
        oszicar = open(str(self.oszicar), "r")
        oszicar_data = oszicar.readlines()

        E0_list = []
        for entry in oszicar_data:
            if 'E0=' in entry:
                E0_list.append(entry)

        E0_datalist = []
        for entry in E0_list:
            E0_data = re.findall("-?0?.?\ *[0-9]+\.?\+?[0-9]*(?:[Ee][E-]\ *-?\ *[0-9]+)?(?:[Ee][E+]\ *-?\ *[0-9]+)?",entry)
            E0_datalist.append(E0_data)

        energy = float(E0_datalist[-1][3])
        oszicar.close()
        return energy

    def get_energy_per_atom(self):
        energy = self.get_energy
        total_atoms = PoscarAnalyzer(self.poscar).get_total_atoms
        energy_per_atom = energy/total_atoms
        return energy_per_atom

    def get_total_magnetization(self):
        oszicar = open(str(self.oszicar), "r")
        oszicar_data = oszicar.readlines()
        mag_list = []
        for entry in oszicar_data:
            if 'mag=' in entry:
                mag_list.append(entry)

        magdata_list = []
        for entry in mag_list:
            mag_data = re.findall(r"[-+]?\d*\.\d+|\d+", entry)
            magdata_list.append(mag_data)

        total_mag = float(magdata_list[-1][-1])
        oszicar.close()
        return total_mag

class DirectoryUtilities(object):

    @staticmethod
    def get_full_directory_list():
        #Get all directories, excluding the cwd
        cwd = os.getcwd()
        directory_list = []
        for directory, subdirs, files in os.walk(cwd):
            if directory != cwd:
                directory_list.append(directory)
        directory_list.sort()
        return directory_list

    @staticmethod
    def get_downmost_directory_list():
        #Get all directories, returning only downmost directories in directory tree
        cwd = os.getcwd()
        directory_list = []
        for directory, subdirs, files in os.walk(cwd):
            if not subdirs:
                directory_list.append(directory)
        directory_list.sort()
        return directory_list

    @staticmethod
    def _parse_directory_list(directory_list_to_parse):
        parsed_directory_list = []
        for entry in directory_list_to_parse:
            # For each directory, check if VASP input files exist. If they do, it's a relevant directory to analyze.
            if os.path.exists(entry+"/"+"POSCAR") and os.path.exists(entry+"/"+"POTCAR") and os.path.exists(entry+"/"+"KPOINTS") and os.path.exists(entry+"/"+"INCAR"):
                parsed_directory_list.append(entry)
            else:
                continue

        return parsed_directory_list

    @staticmethod
    def _modify_directory_list(directory_list, *args):
        # Remove directories from directory_lists based on directory names provided in *args
        modified_directory_list = []; directories_to_remove = []
        for directory_set in args:
            for directory in directory_set:
                directories_to_remove.append(directory)
        for directory in directory_list:
            if directory not in directories_to_remove:
                modified_directory_list.append(directory)
        return modified_directory_list

class TimeUtilities(object):

    @staticmethod
    def _parse_updated_time(last_updated_time):
        last_updated_month = last_updated_time[4:7]
        last_updated_day = int(last_updated_time[8:10])
        last_updated_hour = int(last_updated_time[11:13])
        last_updated_minute = int(last_updated_time[14:16])
        last_updated_second = int(last_updated_time[17:19])
        last_updated_year = int(last_updated_time[20:24])
        if str(last_updated_month) == "Jan":
            last_updated_month = int(1)
        elif str(last_updated_month) == "Feb":
            last_updated_month = int(2)
        elif str(last_updated_month) == "Mar":
            last_updated_month = int(3)
        elif str(last_updated_month) == "Apr":
            last_updated_month = int(4)
        elif str(last_updated_month) == "May":
            last_updated_month = int(5)
        elif str(last_updated_month) == "Jun":
            last_updated_month = int(6)
        elif str(last_updated_month) == "Jul":
            last_updated_month = int(7)
        elif str(last_updated_month) == "Aug":
            last_updated_month = int(8)
        elif str(last_updated_month) == "Sep":
            last_updated_month = int(9)
        elif str(last_updated_month) == "Oct":
            last_updated_month = int(10)
        elif str(last_updated_month) == "Nov":
            last_updated_month = int(11)
        elif str(last_updated_month) == "Dec":
            last_updated_month = int(12)
        last_updated_time_parsed = []
        last_updated_time_parsed.append(last_updated_year); last_updated_time_parsed.append(last_updated_month)
        last_updated_time_parsed.append(last_updated_day); last_updated_time_parsed.append(last_updated_hour)
        last_updated_time_parsed.append(last_updated_minute); last_updated_time_parsed.append(last_updated_second)
        return last_updated_time_parsed

    @staticmethod
    def _get_current_time_data():
        current_time_data = []
        current_time_year = time.strftime("%Y"); current_time_month = time.strftime("%m"); current_time_day = time.strftime("%d")
        current_time_hour = time.strftime("%H"); current_time_minute = time.strftime("%M"); current_time_second = time.strftime("%S")
        current_time_data.append(int(current_time_year)); current_time_data.append(int(current_time_month)); current_time_data.append(int(current_time_day))
        current_time_data.append(int(current_time_hour)); current_time_data.append(int(current_time_minute)); current_time_data.append(int(current_time_second))
        return current_time_data

    @staticmethod
    def _get_file_update_time(directory, filename):
        file_update_time = str(time.ctime(os.path.getmtime(directory + "/" + filename)))
        current_time_data = TimeUtilities._get_current_time_data()
        file_update_time_parsed = TimeUtilities._parse_updated_time(last_updated_time=file_update_time)
        file_update_time_difference = abs(file_update_time_parsed[2] - current_time_data[2]) * 24 + \
                                 abs(file_update_time_parsed[3] - current_time_data[3])
        return file_update_time, file_update_time_difference

class JobAnalyzer(object):

    def __init__(self):
        pass

    def _get_running_and_queued_jobs(self):
        parent_directory = os.getcwd()
        #directory_list = DirectoryUtilities.get_downmost_directory_list()
        running_jobs = []
        queued_jobs = []

        # Get current job log data
        user_id = getpass.getuser()
        job_numbers = []
        try:
            job_status = subprocess.Popen(['squeue', '-u', '%s' % (user_id)], stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                                              stderr=subprocess.PIPE)
        except(OSError):
            # If squeue fails, try qstat.
            job_status = subprocess.Popen(['qstat', '-u', '%s' % (user_id)], stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                                              stderr=subprocess.PIPE)

        # Get the directories that have jobs running, based on their job numbers
        count = 0
        line_number = 0
        jobs_are_running = True
        for line in job_status.stdout.readlines():
            if line_number > 0:
                if "---" in line:
                    continue
                else:
                    job_numbers.append(int(line.split()[0]))
                    count += 1
            line_number += 1
        if count == 0:
            jobs_are_running = False

        if jobs_are_running == False:
            running_and_queued_job_dirs = []
        elif jobs_are_running == True:
            # Find directories which have those job numbers
            for job_number in job_numbers:
                path_with_job_running = subprocess.Popen(['find', parent_directory, '-name', '*%i*' % (job_number)], stdout=subprocess.PIPE)
                if path_with_job_running is not '':
                    directory_with_job_running = os.path.dirname(path_with_job_running.stdout.readline())
                    running_jobs.append(directory_with_job_running)
                else:
                    directory_with_job_running = os.path.dirname(path_with_job_running.stdout.readline())
                    queued_jobs.append(directory_with_job_running)
        """
        # Differentiate between running jobs and queued jobs based on if file containing jobid is written in dir
        running_job_dirs = []
        queued_job_dirs = []

        if len(job_numbers) > 0:
            for directory in running_and_queued_job_dirs:
                for job_number in job_numbers:
                    if job_number in os.listdir(directory):
                        # Job is running
                        running_job_dirs.append(directory)
                    else:
                        # Job is still queued
                        queued_job_dirs.append(directory)

        running_job_dirs_inpath = []
        queued_job_dirs_inpath = []
        for dir1, dir2 in zip(running_job_dirs, queued_job_dirs):
            if dir1 in directory_list:
                running_job_dirs_inpath.append(dir1)
            if dir2 in directory_list:
                queued_job_dirs_inpath.append(dir2)
        """

        #return running_job_dirs_inpath, queued_job_dirs_inpath
        return running_jobs, queued_jobs

    def _get_complete_and_incomplete_jobs(self, directory_list, dE_tolerance=float(10**-2)):
        completed_job_dirs = []; incomplete_job_dirs = []
        #directory_list = DirectoryUtilities.get_downmost_directory_list()

        for directory in directory_list:
            os.chdir(directory)
            try:
                oszicar = open("OSZICAR", "r")
                oszicar_exists = True
            except(IOError):
                continue

            if oszicar_exists == bool(True):
                dE_list = []
                for line in oszicar:
                    if 'd E =' in line:
                        dE_list.append(line)
                dE_datalist = []
                for entry in dE_list:
                    dE_data = re.findall("=?0?-?0?.?\ *[0-9]+\.?\+?[0-9]*(?:[Ee][E-]\ *-?\ *[0-9]+)?(?:[Ee][E+]\ *-?\ *[0-9]+)?", entry)
                    dE_datalist.append(dE_data)
                try:
                    if '=' in dE_datalist[-1][4]:
                        dE_datalist[-1][4] = re.sub('=', '', dE_datalist[-1][4])
                    dE = float(dE_datalist[-1][4])

                    if abs(dE) < abs(dE_tolerance):
                        # In addition to converged energy, check that DOSCAR is written based on file size
                        if os.path.getsize(directory+"/"+"DOSCAR") > 1000:
                            completed_job_dirs.append(directory)
                    elif abs(dE) >= abs(dE_tolerance):
                        # Job most likely finished if the OUTCAR was updated longer than outcar_update_time_tolerance ago
                        outcar_time, outcar_time_difference = TimeUtilities._get_file_update_time(directory=directory, filename="OUTCAR")
                        if outcar_time_difference > 2:
                            # Don't include directories as "incomplete" if they have continuation directories present.
                            if os.path.exists(directory + "/" + "continue"):
                                continue
                            else:
                                incomplete_job_dirs.append(directory)
                        continue

                except(IndexError):
                    continue
                oszicar.close()

        return completed_job_dirs, incomplete_job_dirs

    def _get_crashed_jobs(self, directory_list):
        crashed_job_dirs = []
        #directory_list = DirectoryUtilities.get_downmost_directory_list()
        for directory in directory_list:
            os.chdir(directory)
            # Check if the OSZICAR file exists
            try:
                oszicar = open("OSZICAR", "r")
                oszicar_exists = True
            except(IOError):
                continue

            if oszicar_exists == True:
                dE_list = []
                for line in oszicar:
                    if 'd E =' in line:
                        dE_list.append(line)

                dE_datalist = []
                for entry in dE_list:
                    dE_data = re.findall("=?0?-?0?.?\ *[0-9]+\.?\+?[0-9]*(?:[Ee][E-]\ *-?\ *[0-9]+)?(?:[Ee][E+]\ *-?\ *[0-9]+)?", entry)
                    dE_datalist.append(dE_data)

                try:
                    if '=' in dE_datalist[-1][4]:
                        dE_datalist[-1][4] = re.sub('=', '', dE_datalist[-1][4])
                except(IndexError):
                    crashed_job_dirs.append(directory)
                    continue

                try:
                    dE = dE_datalist[-1][4]
                except(IOError):
                    crashed_job_dirs.append(directory)
                    continue

        return crashed_job_dirs

    def _get_nonstarted_and_old_jobs(self, directory_list):
        directory_time_last_updated_dict = {}
        nonstarted_job_dirs = []
        old_job_dirs = []
        current_job_dirs = []
        for directory in directory_list:
            # Check if the OSZICAR file exists. If it does, get time it was last updated.
            if os.path.exists(directory + "/" + "OSZICAR"):
                last_updated_time = str(time.ctime(os.path.getmtime(directory + "/" + "OSZICAR")))
                last_updated_time_parsed = TimeUtilities._parse_updated_time(last_updated_time=last_updated_time)
                directory_time_last_updated_dict[directory] = last_updated_time_parsed
            # If OSZICAR file doesn't exist, job hasn't run yet.
            else:
                nonstarted_job_dirs.append(directory)

        # Get current time
        current_time_data = TimeUtilities._get_current_time_data()
        update_time_cutoff = 672
        # Check if last update times in directory dictionary are within certain timespan (e.g., 24 hours)
        for key, value in directory_time_last_updated_dict.items():
            # Eliminate directories that were not updated within the update_time_cutoff cutoff
            # Eliminating directories updated in different year
            if value[0] != current_time_data[0]:
                # Eliminate directories updated more than update_time_cutoff ago
                if abs(value[2] - current_time_data[2]) > (update_time_cutoff / 24):
                    old_job_dirs.append(key)
            # If directory updated within value specified by update_time_cutoff, add to current directory list
            # Ensuring directory updated on same year, then check day cutoff, then hour cutoff
            if value[0] == current_time_data[0]:
                if abs(value[2] - current_time_data[2]) <= (update_time_cutoff / 24):
                    if abs((value[3] - current_time_data[3]) * 24) <= update_time_cutoff:
                        current_job_dirs.append(key)

        return nonstarted_job_dirs, old_job_dirs, current_job_dirs

class JobMonitor(JobAnalyzer, DirectoryUtilities, TimeUtilities):

    def __init__(self):
        super(JobMonitor, self).__init__()

    def _write_job_status_report(self, parent_directory, crashed_job_dirs, running_job_dirs, incomplete_job_dirs,
                                completed_job_dirs, resubmitted_job_dirs, submitted_new_job_dirs, nonstarted_job_dirs,
                                old_job_dirs):

        os.chdir(parent_directory)
        job_status_report_file = open("job_status_report.txt", "w")

        # Write directories of crashed jobs to separate file so they can easily be found, only if crashed jobs exist
        if len(crashed_job_dirs) > 0:
            crashed_jobs_file = open("jobs_that_crashed.txt", "w")
            crashed_jobs_file.write("The following jobs crashed and require further attention" + "\n")
            crashed_jobs_file.write("\n")
            for directory in crashed_job_dirs:
                crashed_jobs_file.write(str(directory)+"\n")
            crashed_jobs_file.close()

        # Write directories of newly submitted jobs to separate file so they can be easily found, only if newly submitted jobs exist
        if len(submitted_new_job_dirs) > 0:
            submitted_jobs_file = open("jobs_newly_submitted.txt", "w")
            submitted_jobs_file.write("The following jobs were previously not started and have now been submitted" + "\n")
            submitted_jobs_file.write("\n")
            for entry in submitted_new_job_dirs:
                submitted_jobs_file.write(str(entry)+"\n")
            submitted_jobs_file.close()

        # Write all job data to full status report
        current_time_data = self._get_current_time_data()
        current_time_str = str(current_time_data[0])+"/"+str(current_time_data[1])+"/"+str(current_time_data[2])+", "+"and "+str(current_time_data[3])+" hours, " + str(current_time_data[4])+" minutes, "+str(current_time_data[5])+" seconds"
        job_status_report_file.write("The following job status report file was generated on %s" % (current_time_str) + "\n")
        job_status_report_file.write("\n")

        job_status_report_file.write("The following directories contain jobs that crashed:" + "\n")
        for entry in crashed_job_dirs:
            job_status_report_file.write(str(entry)+"\n")
        if len(crashed_job_dirs)==0:
            job_status_report_file.write("(There are no jobs that were found to have crashed)" + "\n")

        job_status_report_file.write("\n")
        job_status_report_file.write("The following directories contain jobs that are still running:" + "\n")
        for entry in running_job_dirs:
            job_status_report_file.write(str(entry)+"\n")
        if len(running_job_dirs)==0:
            job_status_report_file.write("(There are no jobs that are currently running)" + "\n")

        job_status_report_file.write("\n")
        job_status_report_file.write("The following directories contain jobs that are incomplete:" + "\n")
        for entry in incomplete_job_dirs:
            job_status_report_file.write(str(entry) + "\n")
        if len(incomplete_job_dirs)==0:
            job_status_report_file.write("(There are no jobs that were found to be incomplete)" + "\n")

        job_status_report_file.write("\n")
        job_status_report_file.write("The following directories contain jobs that completed:" + "\n")
        for entry in completed_job_dirs:
            job_status_report_file.write(str(entry)+"\n")
        if len(completed_job_dirs)==0:
            job_status_report_file.write("(There are no jobs that were found to be completed)" + "\n")

        job_status_report_file.write("\n")
        job_status_report_file.write("The following jobs were resubmitted under these directories:" + "\n")
        for entry in resubmitted_job_dirs:
            job_status_report_file.write(str(entry) + "\n")
        if len(resubmitted_job_dirs)==0:
            job_status_report_file.write("(There are no jobs to be resubmitted)" + "\n")

        job_status_report_file.write("\n")
        job_status_report_file.write("The following jobs in these directories have not yet started:" + "\n")
        for entry in nonstarted_job_dirs:
            job_status_report_file.write(str(entry) + "\n")
        if len(nonstarted_job_dirs)==0:
            job_status_report_file.write("(There are no jobs that have not yet started)" + "\n")

        job_status_report_file.write("\n")
        job_status_report_file.write("The following directories contain old jobs that were not analyzed further:" + "\n")
        for entry in old_job_dirs:
            job_status_report_file.write(str(entry) + "\n")
        if len(old_job_dirs)==0:
            job_status_report_file.write("(There are no directories containing old jobs)" + "\n")

        job_status_report_file.close()
        return None

    @staticmethod
    def _resubmit_incomplete_jobs(incomplete_job_dirs):
        resubmitted_job_dirs = []
        for directory in incomplete_job_dirs:
            os.chdir(directory)
            continue_dir = directory+"/"+"continue"
            try:
                os.mkdir(continue_dir)
                shutil.copy("CONTCAR", continue_dir)
                shutil.copy("KPOINTS", continue_dir)
                shutil.copy("INCAR", continue_dir)
                shutil.copy("POTCAR", continue_dir)
                shutil.copy("submit.sh", continue_dir)
                os.chdir(continue_dir)
                shutil.move("CONTCAR", "POSCAR")
                logging.info("Making new continuation folder and resubmitting the job from directory %s" % (continue_dir))
                subprocess.Popen(['sbatch', 'submit.sh']).communicate()
                resubmitted_job_dirs.append(continue_dir)
            except(OSError):
                continue
        return resubmitted_job_dirs

    @staticmethod
    def _submit_nonstarted_jobs(nonstarted_job_dirs):
        newsubmitted_job_dirs = []
        for directory in nonstarted_job_dirs:
            os.chdir(directory)
            logging.info("Submitting new job from directory %s" % (directory))
            try:
                subprocess.Popen(['sbatch', 'submit.sh']).communicate()
            except(OSError):
                subprocess.Popen(['qsub', 'submit.sh']).communicate()
            newsubmitted_job_dirs.append(directory)
        return newsubmitted_job_dirs

    @staticmethod
    def _resubmit_crashed_jobs(crashed_job_dirs):
        resubmitted_crashed_job_dirs = []
        for directory in crashed_job_dirs:
            os.chdir(directory)
            continue_dir = directory+"/"+"continue"
            try:
                os.mkdir(continue_dir)
                shutil.copy("POSCAR", continue_dir)
                shutil.copy("KPOINTS", continue_dir)
                shutil.copy("INCAR", continue_dir)
                shutil.copy("POTCAR", continue_dir)
                shutil.copy("submit.sh", continue_dir)
                os.chdir(continue_dir)
                logging.info("Making new continuation folder and resubmitting the job from directory %s" % (continue_dir))

                # Incar modification for crashed jobs

                incartemp = open("INCARtemp", "w")
                incar = open("INCAR", "r")
                for line in incar:
                    if "NELMDL" not in line and "ALGO" not in line:
                        incartemp.write(line)
                incartemp.write("ALGO = Damped"+"\n")
                #incartemp.write("TIME = 0.05"+"\n")
                incartemp.close()
                incar.close()
                shutil.move("INCARtemp", "INCAR")
                subprocess.Popen(['sbatch', 'submit.sh']).communicate()
                resubmitted_crashed_job_dirs.append(continue_dir)

            except(OSError):
                continue
        return resubmitted_crashed_job_dirs

class VASPDataCollector(object):
    """This class is used to collecting basic VASP output data into a spreadsheet.
    args:
        poscar: (str), the name of a POSCAR file
        outcar: (str), the name of an OUTCAR file
        oszicar: (str), the name of an OSZICAR file
        parent_directory: (str), the name of the parent directory from which to conduct the analysis. If None, will use
        current working directory
    Example usage:
        dc = VASPDataCollector(poscar="POSCAR", outcar="OUTCAR", oszicar="OSZICAR", parent_directory="/home/user/myjobs/")
        dc.get_data_spreadsheet
    """
    def __init__(self, poscar="POSCAR", outcar="OUTCAR", oszicar="OSZICAR", parent_directory=None):
        self.poscar = poscar
        self.outcar = outcar
        self.oszicar = oszicar
        self.parent_directory = parent_directory
        if parent_directory is None:
            parent_directory = os.getcwd()

    @property
    def get_data_spreadsheet(self):
        """
        <Returns>: Instead of returning an object, this method creates a spreadsheet in the working directory containing
        the directory name, composition info, and basic VASP info such as total atoms, energy, fermi energy, dipole,
        cell volume, etc. and saves the file in the directory it is run in.
        <Type>: None, writes xlsx file
        """
        return self._collect_job_data_to_spreadsheet()

    def _collect_job_data_to_spreadsheet(self):
        os.chdir(self.parent_directory)
        excel_file = xlsxwriter.Workbook('job_data_collected.xlsx')
        excel_sheet = excel_file.add_worksheet('collected data')
        bold = excel_file.add_format({'bold': True})
        dirlist = DirectoryUtilities()._get_downmost_directory_list()
        column_name_dict = {"0" : "Run", "1" : "Material", "2" : "Number of atoms", "3" : "Energy (eV/cell)",
                            "4" : "Fermi level (eV)", "5" : "Dipole (eV-A)", "6" : "Volume (A^3)"}
        row = 0
        for key, value in column_name_dict.iteritems():
            #print row, key, value
            excel_sheet.write(int(row), int(key), value, bold)

        row = 1
        for directory in dirlist:
            #if directory != self.parent_directory:
            #print directory
            #print row
            os.chdir(directory)

            # Detect if a directory has the needed files. If not, just pass through:
            try:
                ps = PoscarAnalyzer(poscar=self.poscar)
                oa = OutcarAnalyzer(outcar=self.outcar)
                oza = OszicarAnalyzer(oszicar=self.oszicar, poscar=self.poscar)
                atom_amounts = ps.get_atom_amounts
                element_names = ps.get_element_names
                total_atoms = ps.get_total_atoms
            except(IOError, IndexError):
                atom_amounts = "n/a"
                element_names = "n/a"
                total_atoms = "n/a"

            run = os.getcwd()
            material = str()

            # Detect if a job successfully wrote the following quantities (may not be written if job crashed):
            try:
                dipole = oa.get_dipole
                fermi = oa.get_fermi_energy
                energy = oza.get_energy
                volume = ps.get_cell_volume
            except(IndexError, IOError):
                dipole = "n/a"
                fermi = "n/a"
                energy = "n/a"
                volume = "n/a"

            for index in range(len(element_names)):
                material += str(element_names[index])+str(atom_amounts[index])
            values_dict = {"0" : run, "1" : material, "2" : total_atoms, "3" : energy,
                "4" : fermi, "5" : dipole, "6" : volume}
            for key, value in values_dict.iteritems():
                #print "writing row, column, value to spreadsheet:", row, key, value
                os.chdir(self.parent_directory)
                excel_sheet.write(int(row), int(key), value)
            row += 1

        excel_file.close()