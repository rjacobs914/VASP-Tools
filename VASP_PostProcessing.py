__author__ = 'Ryan Jacobs'
__version__ = '2.0'
__date__ = 'Last updated March 27, 2017'
"""
VASP_PostProcessing is a module designed to handle an assortment of useful post-processing tasks and calculations
based on completed DFT data, including: DOS analysis, thermodynamic stability, Wulff construction, calculation of
diffusion coefficients, and more.
"""

import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot
from pymatgen.analysis import wulff
from pymatgen.core import lattice
from pymatgen.analysis.diffusion_analyzer import DiffusionAnalyzer
from pymatgen.io.vasp.outputs import Xdatcar
from pymatgen.matproj.rest import MPRester
from pymatgen.phasediagram.entries import PDEntry
from pymatgen.core.composition import Composition
from pymatgen.phasediagram.maker import PhaseDiagram
from pymatgen.phasediagram.analyzer import PDAnalyzer
from VASP_Analyzer import PoscarAnalyzer, OutcarAnalyzer, OszicarAnalyzer
import warnings
import math
import numpy as np
import re
from scipy import integrate
import xlsxwriter

class WulffAnalyzer(object):
    """Class used to conduct Wulff construction analysis. Mainly relies on Wulff tools from pymatgen
    args:
        poscar: (str) The name of a POSCAR file
        miller_indices: (list of tuples) A list of miller indices, where each miller index is represented as a tuple.
            For example, miller_indices = [(1, 0, 0), (1, 1, 0), (1, 1, 1)]
        surface_energies: (list of floats) A list of calculated surface energies to pair with the list of miller indices.
    instance methods:
        get_wulff_construction : (wulff object) Creates a wulff construction and saves an image of it to a file
        get_wulff_area_fraction : (dict) Returns the area fraction of each miller index in the Wulff construction
        get_wulff_surface_energies : (dict) Returns the surface energy of each miller index in the Wulff construction
    """
    def __init__(self, miller_indices, surface_energies, poscar="POSCAR"):
        self.miller_indices = miller_indices
        self.surface_energies = surface_energies
        self.poscar = poscar

    def get_wulff_construction(self):
        #Need to get lattice object from POSCAR file
        lattice_parameters = PoscarAnalyzer(poscar=self.poscar).get_lattice_parameters
        lattice_object = lattice.Lattice(lattice_parameters)
        #Make Wulff construction
        wulff_construction = wulff.WulffShape(lattice=lattice_object, miller_list=self.miller_indices,
                                              e_surf_list=self.surface_energies)
        wulff_plot = wulff_construction.get_plot()
        wulff_plot.savefig('wulff_construction.pdf')
        return wulff_construction

    def get_wulff_area_fraction(self):
        wulff_construction = self.get_wulff_construction
        wulff_area = wulff_construction.area_fraction_dict
        return wulff_area

    def get_wulff_surface_energies(self):
        wulff_construction = self.get_wulff_construction
        wulff_energies = wulff_construction.miller_energy_dict
        return wulff_energies

class LocpotAnalyzer(object):
    """This class is designed plot the LOCPOT file, which contains data of the electrostatic potential as a function of
    position in the simulated supercell.
    args:
        poscar: (str) The name of a POSCAR file
        outcar: (str) The name of an OUTCAR file
        locpot: (str) the name of a LOCPOT file
    instance methods:
        get_electrostatic_potential : (list) Returns a list of the planar averaged electrostatic potential from LOCPOT
            file, and also saves the data to an excel file and saves a plot of the data.
        get_workfunction : (tuple of floats) Returns a tuple containing the top surface and bottom surface work functions
        get_vacuum_energy : (float) Returns the vacuum energy of either the top or bottom surface.
            args:
                surface : (str) Set as "Top" to calculate vacuum energy of top surface, analogously for "Bottom"
        get_empirical_delta_workfunction : (float) Returns the difference in top and bottom surface work functions as
            calculated using the VASP calculated dipole and the Helmholtz equation from electrostatics.
    """
    def __init__(self, poscar="POSCAR", outcar="OUTCAR", locpot="LOCPOT"):
        self.poscar = poscar
        self.outcar = outcar
        self.locpot = locpot

    def get_electrostatic_potential(self):
        #Plot the electrostatic potential and save to a file
        pyplot.figure(figsize=(8,6), dpi=80)
        Nx, Ny, Nz, z_coord, LOCPOTdata_column = self._parse_locpot()
        planaravg, avg_planaravg = self._calc_planaravg()
        pyplot.plot(z_coord, planaravg)
        pyplot.xlabel('Z-coordinate (arbitrary units)')
        pyplot.ylabel('Planar averaged electrostatic potential (eV)')
        pyplot.title('Electrostatic potential')
        pyplot.savefig('ElectrostaticPotential.pdf')

        #Output the raw ESP data to an excel file
        excel_file = xlsxwriter.Workbook('Electrostatic_potential_data.xlsx')
        excel_sheet = excel_file.add_worksheet('ESP data')
        bold = excel_file.add_format({'bold': True})
        column_name_dict = {"0" : "Supercell z-coordinate (arb units)", "1" : "Electrostatic potential (eV)"}
        row = 0
        for key, value in column_name_dict.iteritems():
            excel_sheet.write(int(row), int(key), value, bold)
        row = 1
        #values_dict = {"0" : z_coord, "1" : planaravg}
        for index in range(len(z_coord)):
            excel_sheet.write(int(row), 0, z_coord[index])
            excel_sheet.write(int(row), 1, planaravg[index])
            row += 1
        return planaravg

    def get_workfunction(self):
        fermi_energy = OutcarAnalyzer(self.outcar).get_fermi_energy()
        vacuum_energy_top = self.get_vacuum_energy(surface="Top")
        vacuum_energy_bottom = self.get_vacuum_energy(surface="Bottom")
        workfunction_top = vacuum_energy_top-fermi_energy
        workfunction_bottom = vacuum_energy_bottom-fermi_energy
        return (workfunction_top, workfunction_bottom, vacuum_energy_top, vacuum_energy_bottom)

    def get_vacuum_energy(self, surface):
        Nx, Ny, Nz, z_coord, LOCPOTdata_column = self._parse_locpot()
        planaravg, avg_planaravg = self._calc_planaravg()
        # These parameters for finding the values of the vacuum level from each surface seem to work well
        # This will determine how flat the region needs to be
        electrostatic_difference_tolerance = 0.05
        # This determines how many data points for the vacuum level will be averaged together
        count_tol = 5
        # This determines how far from the supercell edges to start looking at values and will determine how large of a flat region needs to be present
        Nz_tol = 10

        vacuum_energy_list = []
        count = 0
        # Loop over range of z-coordinate electrostatic potential values
        if surface == "Top":
            for index in range(Nz):
                # Only consider first count_tol values to average electrostatic potential in vacuum region
                if count < count_tol:
                    # Only consider Nz values that are Nz-Nz_tol away from top supercell edge and Nz_tol away from bottom supercell edge.
                    if index < Nz-Nz_tol and index > Nz_tol:
                        # Only consider electrostatic planaravg values that differ by less than electrostatic_difference_tolerance and are positive.
                        if abs(planaravg[index+Nz_tol]-planaravg[index-Nz_tol]) < electrostatic_difference_tolerance and planaravg[index] > 0:
                            vacuum_energy_list.append(planaravg[index])
                            count += 1
        elif surface == "Bottom":
            Nz_reverse = []
            for index in range(Nz):
                Nz_reverse.append(Nz-index)
            for entry in Nz_reverse:
                if count < count_tol:
                    if entry < Nz-Nz_tol and entry > Nz_tol:
                        if abs(planaravg[entry+Nz_tol]-planaravg[entry-Nz_tol]) < electrostatic_difference_tolerance and planaravg[entry] > 0:
                            vacuum_energy_list.append(planaravg[entry])
                            count += 1
        else:
            raise ValueError('"surface" must be set to either "Top" or "Bottom"')

        sum_vacuum_energy = 0
        for entry in vacuum_energy_list:
            sum_vacuum_energy += entry

        # Remove outliers to try to get a more reasonable vacuum level
        vacuum_energy_top = sum_vacuum_energy/len(vacuum_energy_list)
        stdev_vac = np.std(vacuum_energy_list)
        for entry in vacuum_energy_list:
            if abs(vacuum_energy_top - entry) > stdev_vac:
                vacuum_energy_list.remove(entry)

        sum_vacuum_energy = 0
        for entry in vacuum_energy_list:
            sum_vacuum_energy += entry
        vacuum_energy = sum_vacuum_energy/len(vacuum_energy_list)

        return vacuum_energy

    def get_empirical_delta_workfunction(self):
        #Obtain the work function difference between the two surfaces using the Helmholtz equation
        #This code assumes the surface is oriented in the c-direction, so that a x b is the surface area
        dipole = OutcarAnalyzer(self.outcar).get_dipole
        lattice_parameters = PoscarAnalyzer(self.poscar).get_lattice_parameters
        area = lattice_parameters[0][0]*lattice_parameters[1][1]
        delta_workfunction = ((-181)*dipole)/area #units of eV, dipole is in eV-Ang, area, in Ang^2
        return delta_workfunction

    def _parse_locpot(self):
        LOCPOT = open(self.locpot, "r")
        LOCPOTdata = LOCPOT.readlines()

        #Eliminate the header lines of the LOCPOT file
        total_atoms = int(PoscarAnalyzer(self.poscar).get_total_atoms)
        dataline = LOCPOTdata[total_atoms+9]
        Nx = int(dataline.split()[0])
        Ny = int(dataline.split()[1])
        Nz = int(dataline.split()[2])

        z_coord = []
        for index in range(Nz):
            z_coord.append(index)

        count = 0
        while count < total_atoms+10:
            LOCPOTdata.pop(0)
            count += 1

        LOCPOTdata_column = []
        for entry in LOCPOTdata:
            split_entry = entry.split()
            for entry in split_entry:
                entry = float(entry)
                LOCPOTdata_column.append(entry)

        LOCPOT.close()
        return Nx, Ny, Nz, z_coord, LOCPOTdata_column

    def _calc_planaravg(self):
        Nx, Ny, Nz, z_coord, LOCPOTdata_column = self._parse_locpot()
        planaravg = []
        for z_index in range(Nz):
            planardata = [] #Need to reset planardata values for each z value
            sum_planardata = 0
            for y_index in range(Ny):
                for x_index in range(Nx):
                    datapoint = LOCPOTdata_column[Ny*Nx*(z_index)+Nx*(y_index)+x_index]
                    planardata.append(datapoint)

            for index in range(len(planardata)):
                sum_planardata = sum_planardata + planardata[index]

            planaravg.append(sum_planardata/len(planardata))

        avg_planaravg = 0
        sum_planaravg = 0
        for entry in planaravg:
            sum_planaravg = sum_planaravg + entry

        avg_planaravg = sum_planaravg/Nz
        return planaravg, avg_planaravg

class DiffusionAnalyzerAIMD(object):
    """This class wraps to the pymatgen tools to analyze the output of Ab Initio Molecular Dynamics (AIMD) runs in order
    to calculate and plot the mean-squared displacement (MSD) and calculate the diffusivity components and total
    diffusivity of a given element.
    args:
        diffusing_species: (str) The name of the element to calculate the diffusion coefficient of
        temperature: (int) The temperature (in K) of the AIMD simulation
        timestep: (int) The timestep (in femtoseconds) of the AIMD simulation
        xdatcar: (str) The name of an XDATCAR file
        plot_msd: (bool) Whether or not to plot the MSD.

        **It is advised that the following variables be left as their default values, but if you change them be sure to
        do your own testing**

        steps_to_ignore: (int) The number of initial time steps to leave out of the calculation, generally the ones
            prior to when equilibrium is reached. 100 is typically a good value.
        smoothing: (str) Determines whether or not to smooth the MSD. Choose from three modes: "max", "constant", or "None".
            Recommended to be set to "max".
        min_obs: (int) This is only used if smoothing = "max". In this mode, it specifies the minimum number of
            observations to have before data are included in the MSD calculation. The suggested default is 30.
        avg_nsteps: (int) This is only used if smoothing = "constant". This determines how many timesteps are averaged over
            when iterating over each timestep. The suggested default is 1000.
    instance methods:
        get_diffusion_analysis_from_xdatcars : (dict) Computes the diffusivity of the species of interest and provides
            data of diffusivity and conductivity of the species. Also, plots the MSD and saves to file.
    """
    def __init__(self, diffusing_species, temperature, timestep, xdatcar="XDATCAR", steps_to_ignore=100,
      smoothing="max", min_obs=30, avg_nsteps=1000, plot_msd=True, loop_directories=False):
        self.diffusing_species = diffusing_species
        self.temperature = temperature
        self.timestep = timestep
        self.steps_to_ignore = steps_to_ignore
        self.smoothing = smoothing
        self.min_obs = min_obs
        self.avg_nsteps = avg_nsteps
        self.plot_msd = plot_msd
        self.loop_directories = loop_directories
        self.xdatcar = xdatcar

    def get_diffusion_analysis_from_xdatcars(self):
        structure_list = self._prepare_xdatcar_from_dirs()
        diffusion_analysis = DiffusionAnalyzer.from_structures(structures = structure_list, specie = str(self.diffusing_species),
            temperature = float(self.temperature), time_step = int(self.timestep), step_skip = 1, #Don't change this
            smoothed = str(self.smoothing), min_obs = int(self.min_obs), avg_nsteps = int(self.avg_nsteps))
        outputdict = diffusion_analysis.get_summary_dict(include_msd_t=True)

        if self.plot_msd == bool(True):
            diffusion_analysis.export_msdt("MSD.csv")
            dt = outputdict["dt"]
            msd = outputdict["msd"]
            pyplot.figure(figsize=(8,6), dpi=80)
            pyplot.plot(dt, msd)
            pyplot.xlabel('time (fs)')
            pyplot.ylabel('MSD (Ang^2)')
            pyplot.title('MSD vs time')
            pyplot.savefig('MSD.pdf')

        #print "D of %s (cm^2/sec): %3.3e +/- %3.3e" % (self.diffusing_species, outputdict["D"], outputdict["D_sigma"])
        #print "S of %s (mS/cm): %3.3e +/- %3.3e" % (self.diffusing_species, outputdict["S"], outputdict["S_sigma"])
        #print "D components of %s (cm^2/sec):" % (self.diffusing_species), outputdict["D_components"]
        #print "D sigma components of %s (cm^2/sec):" % (self.diffusing_species), outputdict["D_components_sigma"]
        #print "Specie %s, temp %3.3fK, timestep %i fs, skipping %i initial structures, smoothing %s with min_obs of %i and avg_nsteps of %i" % (self.diffusing_species, self.temperature, self.timestep, self.steps_to_ignore, self.smoothing, self.min_obs, self.avg_nsteps)
        return outputdict

    def _prepare_xdatcar_from_dirs(self):
        structure_list = []
        myxdatcar = Xdatcar(self.xdatcar)
        structure_list.extend(myxdatcar.structures)
        count = 0
        print "Removing %i initial structures from structure list." % self.steps_to_ignore
        while count < self.steps_to_ignore:
            structure_list.pop(0)
            count += 1
        print "%i structures left to analyze." % len(structure_list)
        print "Running diffusion analyzer."
        structures = open("structurelist.txt", "w")
        structures.write(str(structure_list))
        structures.close()
        return structure_list

class DoscarAnalyzer(object):
    """This class is used to obtain plots of projected and total densities of states, and also to obtain useful
    electronic structure quantities from the DOS data, such as bandgap and band centers (centroid of projected DOS for
    certain elements).
    args:
        poscar: (str) the name of a POSCAR file
        incar: (str) the name of an INCAR file
        outcar: (str) the name of an OUTCAR file
        doscar: (str) the name of a DOSCAR file
        energy_cutoff: (float) A user-specified cutoff energy to select a special energy for band center calculations.
            Typically you won't use this input flag.
    instance methods:
        plot_total_dos : (None) saves a plot of the total DOS to file
        plot_projected_dos : (None) saves a plot of the projected DOS to file
        get_bandgap_from_dos : (float) gives the bandgap of the system calculated from the DOS
        get_O_chargetransfergap : (float) gives the charge transfer gap between occupied and unoccupied O states. Only
            relevant for systems containing O.
        get_bandcenters : (list of dict) returns a list of dicts containing the band centers of each orbital type (s, p, etc.)
            of each element in the system.
                args:
                    write_dicts_to_file : (bool) whether to write the bandcenter dicts to a text file.
        get_effective_dos : (tuple of floats) returns the effective density of states for the valence and conduction bands,
            as well as the integral of the valence and conduction bands.
                args:
                    temperature : (int) The effective DOS is T-dependent. Specify a T value (in K) to do the analysis.
    """
    def __init__(self, poscar="POSCAR", incar="INCAR", outcar="OUTCAR", doscar="DOSCAR", energy_cutoff=0):
        self.poscar = poscar
        self.incar = incar
        self.outcar = outcar
        self.doscar = doscar
        self.energy_cutoff = energy_cutoff

    def plot_total_dos(self):
        total_dos_up_list, total_dos_down_list, dos_s_list, dos_p_list, dos_d_list = self._parse_dos_atom_type()
        energy, energy_unocc, energy_occ, index_Fermi, index_cutoff = self._make_energy_lists()
        dosdata, total_dos = self._cleanup_projected_dos()
        number_of_atoms = PoscarAnalyzer(poscar=self.poscar).get_total_atoms

        total_dos_down_list_neg = []
        for index in range(len(total_dos_down_list)):
            total_dos_down_list_neg.append(-1*total_dos_down_list[index])

        total_dos_normalized = []
        for index in range(len(total_dos)):
            total_dos_normalized.append(total_dos[index]/number_of_atoms)

        pyplot.figure(figsize=(8,6), dpi=80)
        pyplot.plot(energy, total_dos_up_list, 'b', label='Total DOS spin up')
        pyplot.plot(energy, total_dos_down_list_neg, 'g', label='Total DOS spin down')
        #pyplot.plot(energy, total_dos_normalized, 'k', label='Total DOS')
        pyplot.plot([0,0], [min(total_dos_up_list), max(total_dos_up_list)], 'k-', label='Fermi energy')
        pyplot.xlabel('Energy (E-EFermi, eV)')
        pyplot.ylabel('Total Densities of States (states/eV-atom)')
        pyplot.title('Total Density of States')
        pyplot.legend()
        pyplot.ylim(-1*max(total_dos), max(total_dos))
        pyplot.xlim([-10, 10])
        pyplot.savefig('TotalDOS.pdf')
        return None

    def plot_projected_dos(self):
        sumtotaldosup, sumtotaldosdown, dos_s_list, dos_p_list, dos_d_list = self._parse_dos_atom_type()
        energy, energy_unocc, energy_occ, index_Fermi, index_cutoff = self._make_energy_lists()
        nedos = self._calc_nedos()
        atom_amounts = PoscarAnalyzer(self.poscar).get_atom_amounts()
        element_names = PoscarAnalyzer(self.poscar).get_element_names()

        pyplot.figure(figsize=(8,6), dpi=80)
        #Need to normalize and plot PDOS by element type

        colors = "bgrkybgrkybgrkybgrkybgrky"
        for element_number in range(len(atom_amounts)):
            dos_to_plot_up = []
            dos_to_plot_down = []
            for index in range(nedos):
                dos_to_plot_up.append((dos_s_list[index + nedos*element_number][0] + dos_p_list[index + nedos*element_number][0] + dos_d_list[index + nedos*element_number][0])/atom_amounts[element_number])
                dos_to_plot_down.append(-1*(dos_s_list[index + nedos*element_number][1] + dos_p_list[index + nedos*element_number][1] + dos_d_list[index + nedos*element_number][1])/atom_amounts[element_number])
            pyplot.plot(energy, dos_to_plot_up, color=colors[element_number], label='PDOS for %s' % (element_names[element_number]))
            pyplot.plot(energy, dos_to_plot_down, color=colors[element_number])

        pyplot.xlabel('Energy (E-EFermi, eV)')
        pyplot.ylabel('Projected Densities of States (states/eV-atom)')
        pyplot.title('Projected Density of States')
        pyplot.legend()
        pyplot.xlim([-10, 10])
        pyplot.ylim([-5, 5])
        pyplot.savefig('ProjectedDOS.pdf')
        return None

    def get_bandgap_from_dos(self):
        E_VBM, E_CBM = self._calc_VBM_and_CBM()
        Egap = E_CBM - E_VBM
        if Egap < 0.1:
            Egap = 0

        #print "The bandgap of this material from DOS is %s eV" % (Egap)
        bandgapfile = open("bandgap_fromdos.txt", "w")
        Egap = str(Egap)
        bandgapfile.write(Egap)
        bandgapfile.close()
        return Egap

    def get_O_chargetransfergap(self):
        # Same method as regular bandgap, but now just doing it between occupied and unoccupied O states
        total_dos_up_list, total_dos_down_list, dos_s_list, dos_p_list, dos_d_list = self._parse_dos_atom_type()
        energy, energy_unocc, energy_occ, index_Fermi, index_cutoff = self._make_energy_lists()
        nedos = self._calc_nedos()
        atom_amounts = PoscarAnalyzer(self.poscar).get_atom_amounts()
        element_names = PoscarAnalyzer(self.poscar).get_element_names()
        Etolerance = 0.05
        DOStolerance = 0.001

        # Find which element, if any, is O. This analysis only pertains to O, so throw error if no O found
        if "O" not in element_names:
            raise TypeError("No oxygen was detected in your element list. The charge transfer gap calculation only pertains to O states. Exiting calculation...")

        # Find out which entry is O and perform charge transfer gap calculation
        for element_number in range(len(element_names)):
            if element_names[element_number] == "O":
                # Add the respective p up and p down states based on the element entry
                # print "O is element number", element_number
                dos_O_p_total = []
                for dos_row in range(nedos):
                    dos_O_p_total.append((dos_p_list[dos_row + nedos * element_number][0] +
                                          dos_p_list[dos_row + nedos * element_number][1]) / atom_amounts[
                                             element_number])

                # Find the position of the VBM
                for index in range(nedos):
                    if energy[index] >= 0 - Etolerance:
                        if energy[index] <= 0 + Etolerance:
                            E_VBM = energy[index]
                            VBM_index = index

                    elif energy[index] < 0 - Etolerance or energy[index] > 0 + Etolerance:
                        pass

                E_above_VBM_range = []
                for index in range(nedos - VBM_index):
                    E_above_VBM_range.append(VBM_index + index)

                # Find the position of the last O state after Efermi (where O states go to zero)
                for entry in E_above_VBM_range:
                    if abs(dos_O_p_total[entry] - dos_O_p_total[entry - 1]) > 0 + DOStolerance:
                        pass
                    elif abs(dos_O_p_total[entry] - dos_O_p_total[entry - 1]) <= 0 + DOStolerance:
                        new_E_VBM = energy[entry]
                        new_VBM_index = entry
                        break

                new_E_above_VBM_range = []
                for index in range(nedos - new_VBM_index):
                    new_E_above_VBM_range.append(new_VBM_index + index)

                for entry in new_E_above_VBM_range:
                    if abs(dos_O_p_total[entry] - dos_O_p_total[entry - 1]) > 0 + DOStolerance:
                        new_CBM_index = entry
                        new_E_CBM = energy[entry]
                        break
                    elif abs(dos_O_p_total[entry] - dos_O_p_total[entry - 1]) < 0 + DOStolerance:
                        pass

                Echgtransgap = new_E_CBM - new_E_VBM
                if Echgtransgap < 0.1:
                    Echgtransgap = 0

                #print "The charge transfer gap for O in this material is %s eV" % (Echgtransgap)
                bandgapfile = open("chgtransgap.txt", "w")
                Echgtransgap = str(Echgtransgap)
                bandgapfile.write(Echgtransgap)
                bandgapfile.close()
                return Echgtransgap

    def get_bandcenters(self, write_dicts_to_file=True):
        nedos = self._calc_nedos()
        energy, energy_unocc, energy_occ, index_Fermi, index_cutoff = self._make_energy_lists()
        sumtotaldosup, sumtotaldosdown, dos_s_list, dos_p_list, dos_d_list = self._parse_dos_atom_type()
        pa = PoscarAnalyzer(self.poscar)
        element_names = pa.get_element_names()
        atom_amounts = pa.get_atom_amounts()
        dos_s_list_E = []
        dos_s_list_total_occ = []
        dos_p_list_E = []
        dos_p_list_total_occ = []
        dos_d_list_E = []
        dos_d_list_total_occ = []
        dos_s_list_total = []
        dos_s_list_E_occ = []
        dos_p_list_total = []
        dos_p_list_E_occ = []
        dos_d_list_total = []
        dos_d_list_E_occ = []

        band_center_list = []
        band_center_occ_list = []

        for element_number in range(len(atom_amounts)):
            for index in range(nedos):
                dos_s_list_total.append(
                    dos_s_list[index + nedos * element_number][0] + dos_s_list[index + nedos * element_number][1])
                dos_p_list_total.append(
                    dos_p_list[index + nedos * element_number][0] + dos_p_list[index + nedos * element_number][1])
                dos_d_list_total.append(
                    dos_d_list[index + nedos * element_number][0] + dos_d_list[index + nedos * element_number][1])
                dos_s_list_E.append(dos_s_list_total[index] * energy[index])
                dos_p_list_E.append(dos_p_list_total[index] * energy[index])
                dos_d_list_E.append(dos_d_list_total[index] * energy[index])
                if index < index_Fermi:
                    dos_s_list_total_occ.append(
                        dos_s_list[index + nedos * element_number][0] + dos_s_list[index + nedos * element_number][1])
                    dos_p_list_total_occ.append(
                        dos_p_list[index + nedos * element_number][0] + dos_p_list[index + nedos * element_number][1])
                    dos_d_list_total_occ.append(
                        dos_d_list[index + nedos * element_number][0] + dos_d_list[index + nedos * element_number][1])
                    dos_s_list_E_occ.append(dos_s_list_total[index] * energy[index])
                    dos_p_list_E_occ.append(dos_p_list_total[index] * energy[index])
                    dos_d_list_E_occ.append(dos_d_list_total[index] * energy[index])

            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                s_bandcenter = integrate.cumtrapz(dos_s_list_E, energy)[-1] / \
                               integrate.cumtrapz(dos_s_list_total, energy)[-1]
                p_bandcenter = integrate.cumtrapz(dos_p_list_E, energy)[-1] / \
                               integrate.cumtrapz(dos_p_list_total, energy)[-1]
                d_bandcenter = integrate.cumtrapz(dos_d_list_E, energy)[-1] / \
                               integrate.cumtrapz(dos_d_list_total, energy)[-1]
                s_bandcenter_occ = integrate.cumtrapz(dos_s_list_E_occ, energy_occ)[-1] / \
                                   integrate.cumtrapz(dos_s_list_total_occ, energy_occ)[-1]
                p_bandcenter_occ = integrate.cumtrapz(dos_p_list_E_occ, energy_occ)[-1] / \
                                   integrate.cumtrapz(dos_p_list_total_occ, energy_occ)[-1]
                d_bandcenter_occ = integrate.cumtrapz(dos_d_list_E_occ, energy_occ)[-1] / \
                                   integrate.cumtrapz(dos_d_list_total_occ, energy_occ)[-1]
                band_center_dict = {}
                band_center_dict_occ = {}
                element = element_names[element_number]
                band_center_dict[str(element) + "_s"] = "%3.3f" % s_bandcenter
                band_center_dict[str(element) + "_p"] = "%3.3f" % p_bandcenter
                band_center_dict[str(element) + "_d"] = "%3.3f" % d_bandcenter
                band_center_dict_occ[str(element) + "_s_occ"] = "%3.3f" % s_bandcenter_occ
                band_center_dict_occ[str(element) + "_p_occ"] = "%3.3f" % p_bandcenter_occ
                band_center_dict_occ[str(element) + "_d_occ"] = "%3.3f" % d_bandcenter_occ

                if write_dicts_to_file == bool(True):
                    # Write bandcenter dicts to text files
                    filename = "%s_bandcenters.txt" % element
                    file = open(filename, "w")
                    for key, value in band_center_dict.items():
                        file.write(str(key)+"="+str(value))
                    file.close()
                    filename = "%s_bandcenters_occ.txt" % element
                    file = open(filename, "w")
                    for key, value in band_center_dict_occ.items():
                        file.write(str(key)+"="+str(value))
                    file.close()

                band_center_list.append(band_center_dict)
                band_center_occ_list.append(band_center_dict_occ)

        return band_center_list, band_center_occ_list

    def get_effective_dos(self, temperature=1000):
        boltz = 8.619 * 10 ** -5 #units of eV/K
        E_VBM, E_CBM = self._calc_VBM_and_CBM()
        energy, energy_unocc, energy_occ, index_Fermi, index_cutoff = self._make_energy_lists()
        dosdata, total_dos = self._cleanup_projected_dos()

        # Parse total DOS to VBM and CBM components
        total_dos_VBM = []
        energy_to_integrate_VBM = []
        total_dos_CBM = []
        energy_to_integrate_CBM = []
        for index in range(len(energy)):
            # if energy[index] <= E_VBM and energy[index] > E_VBM-10:
            if energy[index] <= E_VBM:
                total_dos_VBM.append(total_dos[index])
                energy_to_integrate_VBM.append(energy[index])
            # elif energy[index] > E_CBM and energy[index] < E_CBM+10:
            elif energy[index] > E_CBM:
                total_dos_CBM.append(total_dos[index])
                energy_to_integrate_CBM.append(energy[index])
            else:
                continue

        # Calculate the effective DOS of the valence band
        VBM_effective_dos_integrand = []
        for index in range(len(total_dos_VBM)):
            VBM_effective_dos_integrand.append(
                total_dos_VBM[index] * np.exp(-1 * (E_VBM - energy_to_integrate_VBM[index]) / (boltz * temperature)))
        VBM_effective_dos = integrate.cumtrapz(VBM_effective_dos_integrand, energy_to_integrate_VBM)[-1]
        VBM_int = integrate.cumtrapz(total_dos_VBM, energy_to_integrate_VBM)[-1]

        # Calculate the effective DOS of the conduction band
        CBM_effective_dos_integrand = []
        for index in range(len(total_dos_CBM)):
            CBM_effective_dos_integrand.append(
                total_dos_CBM[index] * np.exp(-1 * (energy_to_integrate_CBM[index] - E_CBM) / (boltz * temperature)))
        CBM_effective_dos = integrate.cumtrapz(CBM_effective_dos_integrand, energy_to_integrate_CBM)[-1]
        CBM_int = integrate.cumtrapz(total_dos_CBM, energy_to_integrate_CBM)[-1]

        #print "The effective DOS of the VB is:", VBM_effective_dos
        #print "The effective DOS of the CB is:", CBM_effective_dos
        #print "The integral of the VB DOS is:", VBM_int
        #print "The integral of the CB DOS is:", CBM_int

        return (VBM_effective_dos, CBM_effective_dos, VBM_int, CBM_int)

    def _calc_nedos(self):
        # Obtain the value of NEDOS from INCAR file
        incar = open(self.incar, "r")
        nedos_found = False
        for line in incar:
            if 'NEDOS' in line:
                nedos_found = True
                NEDOS_info = line
                NEDOS_info = re.findall('\d+',NEDOS_info)
                for entry in NEDOS_info:
                    nedos = int(entry)

        #If no NEDOS is specified in INCAR, use default value of 301
        if nedos_found == bool(False):
            nedos = int(301)

        incar.close()
        return nedos

    def _make_energy_lists(self):
        Energy = []
        dosdata, total_dos = self._cleanup_projected_dos()
        nedos = self._calc_nedos()
        fermi = OutcarAnalyzer(self.outcar).get_fermi_energy()
        # Create the list of energies
        count = 0
        for entry in dosdata:
            if count < nedos:
                Energy.append(float(entry[0:12])-fermi)
            count += 1

        count = 0
        for entry in Energy:
            tol = 0.05
            count += 1
            if entry > 0-tol and entry < 0+tol: #So the Energy should be between -0.2 and 0.2 eV
                index_Fermi = count

        count = 0
        if self.energy_cutoff is not None:
            for entry in Energy:
                tol = 0.05
                count += 1
                if entry > self.energy_cutoff-tol and entry < self.energy_cutoff+tol:
                    index_cutoff = count

        Energy_occ = []
        Energy_unocc = []
        for index in range(nedos):
            if index < index_Fermi:
                Energy_occ.append(Energy[index])
            if index > index_Fermi:
                Energy_unocc.append(Energy[index])

        # This is used to create a custom energy cutoff that is not at the Fermi level
        Energy_range_below = []
        Energy_range_above = []
        if self.energy_cutoff is not None:
            for index in range(nedos):
                if index < index_cutoff:
                    Energy_range_below.append(Energy[index])
                if index > index_cutoff:
                    Energy_range_above.append(Energy[index])

        return Energy, Energy_unocc, Energy_occ, index_Fermi, index_cutoff

    def _cleanup_projected_dos(self):
        # Remove the header lines of the DOSCAR file
        doscar = open(self.doscar, "r")
        dosdata = doscar.readlines()
        nedos = self._calc_nedos()
        number_of_atoms = int(PoscarAnalyzer(self.poscar).get_total_atoms())

        # Isolate just the total DOS from the DOS data
        total_dos_rows = []
        total_dos = []
        count = 0
        for line in dosdata:
            if count > 5 and count < nedos+6:
                total_dos_rows.append(line)
            count += 1

        for index in range(len(total_dos_rows)):
            #print float(total_dos_rows[index].split()[1])
            total_dos.append(float(total_dos_rows[index].split()[1]) + float(total_dos_rows[index].split()[2]))

        count = 0
        while count < nedos + 7:
            dosdata.pop(0)
            count += 1

        # Eliminate the interspersed junk lines between each atom in DOS data
        for index in range(len(dosdata)-number_of_atoms):
            if index % nedos == 0 and index != 0:
                dosdata.pop(index)

        doscar.close()
        return dosdata, total_dos

    def _parse_projected_dos_string_to_float(self):
        dosdata_fromfile, total_dos_fromfile = self._cleanup_projected_dos()

        dosdata = []
        #Split the lines of the doscar file and make all entries floats rather than strings
        for line in dosdata_fromfile:
            line = line.split()
            for index in range(len(line)):
                line[index] = float(line[index])
            dosdata.append(line)

        return dosdata

    def _parse_dos_orbital_type(self):
        dosdata = self._parse_projected_dos_string_to_float()
        doscolumns = len(dosdata[0])
        dosdata_orbital_summed = []
        dosdata_orbital_summed_float = []
        for dosrow in range(len(dosdata)):
            sum_s_up = 0
            sum_s_down = 0
            sum_p_up = 0
            sum_p_down = 0
            sum_d_up = 0
            sum_d_down = 0
            for doscolumn in range(doscolumns):
                if doscolumn > 0: #Don't include values of energy in sum, only want to sum DOS values
                    if doscolumn == 1: #s up
                        sum_s_up += dosdata[dosrow][doscolumn]
                    if doscolumn == 2:
                        sum_s_down += dosdata[dosrow][doscolumn]
                    if doscolumn == 3:
                        sum_p_up += dosdata[dosrow][doscolumn] + dosdata[dosrow][doscolumn+2] + dosdata[dosrow][doscolumn+4] #Adding columns 3, 5, and 7, which are px, py and pz up
                    if doscolumn == 4:
                        sum_p_down += (dosdata[dosrow][doscolumn] + dosdata[dosrow][doscolumn+2] + dosdata[dosrow][doscolumn+4]) #Adding columns 4, 6, and 8, which are px, py and pz down
                    if doscolumn > 4 and doscolumn < 9:
                        pass
                    if doscolumn == 9:
                        sum_d_up += dosdata[dosrow][doscolumn] + dosdata[dosrow][doscolumn+2] + dosdata[dosrow][doscolumn+4] + dosdata[dosrow][doscolumn+6] + dosdata[dosrow][doscolumn+8] #Adding columns 9, 11, 13, 15 and 17 which are dxy, dyz, dz2, dxz, dx2 up
                    if doscolumn == 10:
                        sum_d_down += (dosdata[dosrow][doscolumn] + dosdata[dosrow][doscolumn+2] + dosdata[dosrow][doscolumn+4] + dosdata[dosrow][doscolumn+6] + dosdata[dosrow][doscolumn+8]) #Adding columns 10, 12, 14, 16 and 18 which are dxy, dyz, dz2, dxz, dx2 down
                    if doscolumn > 10:
                        pass
            dosdata_orbital_summed.append(str(sum_s_up)+" "+str(sum_s_down)+" "+str(sum_p_up)+" "+str(sum_p_down)+" "+str(sum_d_up)+" "+str(sum_d_down)+"\n")

        for line in dosdata_orbital_summed:
            line = line.split()
            for index in range(len(line)):
                line[index] = float(line[index])
            dosdata_orbital_summed_float.append(line)

        return dosdata_orbital_summed_float

    def _parse_dos_atom_type(self):
        dosdata = self._parse_dos_orbital_type()
        nedos = self._calc_nedos()
        atom_amounts = PoscarAnalyzer(self.poscar).get_atom_amounts()
        element_names = PoscarAnalyzer(self.poscar).get_element_names()

        element_iteration_count = 0
        dos_sum_s_up_list = []
        dos_sum_s_down_list = []
        dos_sum_p_up_list = []
        dos_sum_p_down_list = []
        dos_sum_d_up_list = []
        dos_sum_d_down_list = []
        for index in range(len(atom_amounts)): #Loop over each atom type
            #print "analyzing atom type", element_names[index]
            #print "element", element_names[index], "is element number", element_iteration_count
            for dos_row in range(nedos):
                dos_sum_s_up = 0
                dos_sum_s_down = 0
                dos_sum_p_up = 0
                dos_sum_p_down = 0
                dos_sum_d_up = 0
                dos_sum_d_down = 0
                for atom_number in range(atom_amounts[index]):
                    if index == 0:
                        true_atom_number = atom_number
                    if index == 1:
                        true_atom_number = atom_number + atom_amounts[index-1]
                    if index == 2:
                        true_atom_number = atom_number + atom_amounts[index-1] + atom_amounts[index-2]
                    if index == 3:
                        true_atom_number = atom_number + atom_amounts[index-1] + atom_amounts[index-2] + atom_amounts[index-3]
                    if index == 4:
                        true_atom_number = atom_number + atom_amounts[index-1] + atom_amounts[index-2] + atom_amounts[index-3] + atom_amounts[index-4]
                    if index == 5:
                        true_atom_number = atom_number + atom_amounts[index-1] + atom_amounts[index-2] + atom_amounts[index-3] + atom_amounts[index-4] + atom_amounts[index-5]
                    dos_sum_s_up += dosdata[dos_row + nedos*true_atom_number][0]
                    dos_sum_s_down += dosdata[dos_row + nedos*true_atom_number][1]
                    dos_sum_p_up += dosdata[dos_row + nedos*true_atom_number][2]
                    dos_sum_p_down += dosdata[dos_row + nedos*true_atom_number][3]
                    dos_sum_d_up += dosdata[dos_row + nedos*true_atom_number][4]
                    dos_sum_d_down += dosdata[dos_row + nedos*true_atom_number][5]
                dos_sum_s_up_list.append(dos_sum_s_up)
                dos_sum_s_down_list.append(dos_sum_s_down)
                dos_sum_p_up_list.append(dos_sum_p_up)
                dos_sum_p_down_list.append(dos_sum_p_down)
                dos_sum_d_up_list.append(dos_sum_d_up)
                dos_sum_d_down_list.append(dos_sum_d_down)
            element_iteration_count += 1

        # The dos_sum lists will contain nedos*number_of_elements lines. So, 4 elements with 2000 nedos will be 8000 lines.
        dos_s_list = []
        dos_p_list = []
        dos_d_list = []
        dos_sum_up_list = []
        dos_sum_down_list = []
        for index in range(len(dos_sum_s_up_list)):
            dos_s_list.append(str(dos_sum_s_up_list[index])+" "+str(dos_sum_s_down_list[index])+"\n")
            dos_p_list.append(str(dos_sum_p_up_list[index])+" "+str(dos_sum_p_down_list[index])+"\n")
            dos_d_list.append(str(dos_sum_d_up_list[index])+" "+str(dos_sum_d_down_list[index])+"\n")
            dos_sum_up_list.append(dos_sum_s_up_list[index] + dos_sum_p_up_list[index] + dos_sum_d_up_list[index])
            dos_sum_down_list.append(dos_sum_s_down_list[index] + dos_sum_p_down_list[index] + dos_sum_d_down_list[index])

        dos_s_list_float = []
        dos_p_list_float = []
        dos_d_list_float = []
        for line in dos_s_list:
            line = line.split()
            for index in range(len(line)):
                line[index] = float(line[index])
            dos_s_list_float.append(line)
        for line in dos_p_list:
            line = line.split()
            for index in range(len(line)):
                line[index] = float(line[index])
            dos_p_list_float.append(line)
        for line in dos_d_list:
            line = line.split()
            for index in range(len(line)):
                line[index] = float(line[index])
            dos_d_list_float.append(line)

        #Now get total DOS by summing all projected dos
        total_dos_up_list = []
        total_dos_down_list = []
        for dos_row in range(nedos):
            sum_dos_up = 0
            sum_dos_down = 0
            for index in range(len(atom_amounts)):
                # These are normalized by number of atoms of each element type
                sum_dos_up += (dos_sum_up_list[dos_row + nedos*index])/atom_amounts[index]
                sum_dos_down += (dos_sum_down_list[dos_row + nedos*index])/atom_amounts[index]
                # These are not normalized by atom number
                #sum_dos_up += (dos_sum_up_list[dos_row + nedos*index])
                #sum_dos_down += (dos_sum_down_list[dos_row + nedos*index])

            total_dos_up_list.append(sum_dos_up) #This is normalized by atom amounts of each element present
            total_dos_down_list.append(sum_dos_down)

        #print "length total dos up", len(total_dos_up_list)
        #print "length sum dos up list", len(dos_sum_up_list)
        return total_dos_up_list, total_dos_down_list, dos_s_list_float, dos_p_list_float, dos_d_list_float

    def _calc_VBM_and_CBM(self):
        total_dos_up_list, total_dos_down_list, dos_s_list, dos_p_list, dos_d_list = self._parse_dos_atom_type()
        energy, energy_unocc, energy_occ, index_Fermi, index_cutoff = self._make_energy_lists()
        nedos = self._calc_nedos()
        Etolerance = 0.1
        DOStolerance = 0.01

        #Take the total dos as separate up and down spin and combine to be summed total DOS
        total_dos = []
        for index in range(len(total_dos_up_list)):
            total_dos.append(total_dos_up_list[index]+total_dos_down_list[index])

        #Find the position of the VBM
        for index in range(nedos):
            if energy[index] >= 0 - Etolerance:
                if energy[index] <= 0 + Etolerance:
                    E_VBM = energy[index]
                    VBM_index = index

            elif energy[index] < 0 - Etolerance or energy[index] > 0 + Etolerance:
                pass

        E_above_VBM_range = []
        for index in range(nedos - VBM_index):
            E_above_VBM_range.append(VBM_index + index)

        #Now calculate the band gap
        for entry in E_above_VBM_range:
            if abs(total_dos[entry] - total_dos[entry-1]) > 0 + DOStolerance:
                CBM_index = entry
                E_CBM = energy[entry]
                break
            elif abs(total_dos[entry] - total_dos[entry-1]) < 0 + DOStolerance:
                pass

        return E_VBM, E_CBM

class StabilityAnalyzer(object):
    """Class that conducts a thermodynamic phase stability analysis using the values in the Materials Project database and
    includes as input the user's calculated material. By default, DFT energies are used for all elemental endmembers.
    However, the thermodynamic conditions represented by DFT are often not physically realistic. Because of this, the user
    can supply a custom set of calculated chemical potential values for gaseous endmembers like O and H. These chemical
    potentials can be calculated using the ChemicalPotentialAnalyzer class.

    * A couple notes about usage *
    Since this entire analysis package uses the Materials Project database, it is important that the same XC functional,
    pseudopotentials and GGA+U values are used so that consistent energies are obtained. If you use the IncarFileSetup
    class in the VASP_Setup module your GGA+U values will match that of the Materials Project. Just be sure to use GGA-PBE
    type pseudopotentials in your DFT calculations.

    args:
        mapi_key : (str) Your Materials API key from Materials Project. You can obtain yours here: https://materialsproject.org/open
        poscar : (str) the name of a POSCAR file
        oszicar : (str) the name of an OSZICAR file
        get_data_from_VASP_files : (bool) whether to obtain composition and energy data from VASP runs. This is used to
            calculate the stability of newly calculated DFT structures. If this is set to False, need to specify composition
            and energy of material you want to examine manually with the composition_dict and composition_energy arguments
        additional_elements_to_include : (list) A list of elements to include in phase stability analysis that aren't present
            in your DFT calculation or that weren't manually included in the composition_dict argument
        material_ids_to_remove : (list) A list of MP material id's (e.g. "mp-123") to remove from the phase space. Useful
            to analyze metastable phases by removing certain stable phases
        composition_dict : (list of dict) A list containing dicts of material compositions to analyze the stability of.
            This is a way to manually set the composition of interest for an arbitrary composition
        composition_energy : (list) A list of energies corresponding to each composition dict in composition_dict

    instance methods:
        get_phase_diagram: (phase diagram object) This function generates a phase diagram and calculates the stability
            of the simulated compound
            args:
                use_custom_chem_pots : (bool) If False, the standard DFT endmember energies as tabulated in Materials Project
                    are used. If True, the user must specify their own chemical potentials using the custom_chem_pot_dict input
                custom_chem_pot_dict: (dict) A dictionary of {"element": chem_pot} pairs, where chem_pot is a float
                include_organic_molecules: (bool) If True, a small library of DFT-calculated organic molecules is included
                    in the phase diagram calculation. Note that they can only be included if C is an element in your POSCAR.
                include_organic_molecule_shift: (bool) If True, adds an energy shift to the organic molecules consistent with
                    T = 298 K and P = 10^-6 atm. Note this only applies if include_organic_molecules is set to True and you
                    are interested in phase stability at 298 K. Use with Caution!
    """
    def __init__(self, mapi_key, poscar="POSCAR", oszicar="OSZICAR", get_data_from_VASP_files=False ,
                 additional_elements_to_include=None, material_ids_to_remove=None, composition_dict=None,
                 composition_energy=None):
        self.mapi_key = mapi_key
        self.poscar = poscar
        self.oszicar = oszicar
        self.get_data_from_VASP_files = get_data_from_VASP_files
        self.additional_elements_to_include = additional_elements_to_include
        self.material_ids_to_remove = material_ids_to_remove
        self.composition_dict = composition_dict
        self.composition_energy = composition_energy


    def get_phase_diagram(self, use_custom_chem_pots=False, custom_chem_pot_dict=None, include_organic_molecules=False,
                          include_organic_molecule_shift=False):

        pd_entry_list_chempots = []

        if use_custom_chem_pots == bool(False):
            entries_in_system = self._create_chemical_system_from_MP()
            pd_entry_list = self._create_pdentry()
            for entry in pd_entry_list:
                print "The PDEntry given for this system is:", entry, " eV/cell"

        elif use_custom_chem_pots == bool(True):
            entries_in_system = self._create_chemical_system_from_MP_and_remove_endmembers(species_to_remove=[key for key in custom_chem_pot_dict.keys()])
            pd_entry_list = self._create_pdentry()
            # Also make PDEntries for species you specify chem pots of
            for key, value in custom_chem_pot_dict.items():
                chempot_pdentry = PDEntry(composition=key, energy=value)
                entries_in_system.append(chempot_pdentry)
                pd_entry_list_chempots.append(chempot_pdentry)
            for entry in pd_entry_list:
                print "The PDEntry given for this system is:", entry, " eV/cell"
            for entry in pd_entry_list_chempots:
                print "The PDEntry of new chemical potential is:", entry, " eV"

        # Add all new PDEntry objects from pd_entry_list to the chemical system for phase stability analysis
        if len(pd_entry_list) > 0:
            for entry in pd_entry_list:
                entries_in_system.append(entry)
        if len(pd_entry_list_chempots) > 0:
            for entry in pd_entry_list_chempots:
                entries_in_system.append(entry)

        # Add the organic molecule entries as PDEntries and entries in chemical system for phase diagram
        if include_organic_molecules == bool(True):
            organic_count = 0
            organic_compositions, organic_energies = self._add_organic_molecules(add_roomtemp_gas_shift=include_organic_molecule_shift)
            for index in range(len(organic_compositions)):
                comp = Composition(organic_compositions[index])
                energy = float(organic_energies[index])
                organic_pdentry = PDEntry(comp, energy)
                pd_entry_list.append(organic_pdentry)
                entries_in_system.append(organic_pdentry)
                organic_count += 1
            print "%i organic molecules have been added to the chemical system!" % organic_count

        phasediagram = PhaseDiagram(entries=entries_in_system)
        phasediagram_analyzer = PDAnalyzer(phasediagram)

        eabove_file = open("energy_above_hull.txt", "w")
        eform_file = open("formation_energy.txt", "w")
        decomp_file = open("decomposition_products.txt", "w")
        stable_entries_file = open("stable_entries.txt", "w")

        for entry in pd_entry_list:
            print "Analyzing entry:", entry
            energy_above_hull = phasediagram_analyzer.get_e_above_hull(entry=entry)
            energy_above_hull *= 1000 # in units of meV/atom
            print "The energy above hull (in meV/atom) for this system is:", energy_above_hull
            eabove_file.write(str(energy_above_hull)+"\n")
            e_form = phasediagram.get_form_energy_per_atom(entry=entry)
            #e_form *= 1000 # in units of meV/atom
            print "The formation energy (in eV/atom) for this system is:", e_form
            eform_file.write(str(e_form)+"\n")

        # Write the stable entries to a file
        for entry in phasediagram.stable_entries:
            if entry not in pd_entry_list and entry not in pd_entry_list_chempots:
                stable_entries_file.write(str(entry.composition.reduced_formula)+" "+str(entry.entry_id)+"\n")
            if entry in pd_entry_list or entry in pd_entry_list_chempots:
                stable_entries_file.write(str(entry)+"\n")

        # Write the decomposition products and fractions of each phase for inserted PDEntry to file
        for entry in phasediagram.unstable_entries:
            if entry in pd_entry_list:
                decomp, e_above_hull = phasediagram_analyzer.get_decomp_and_e_above_hull(entry)
                for key, value in decomp.iteritems():
                    decomp_file.write(str(key.composition.reduced_formula)+" "+str(value)+"\n")
                decomp_file.write("\n")

        eabove_file.close()
        eform_file.close()
        decomp_file.close()
        stable_entries_file.close()

        return phasediagram

    def _calc_GGA_GGAU_shifted_energy(self, composition=None):
        # A note of caution for this function and all stability analysis: for your answers to be meaningful, you need
        # to (1) use the same XC-functional as MP (PAW-PBE), and (2) use same U values for GGA+U calculations. The
        # shifts done here correspond to certain U values used by MP, which are the same U values used in the
        # VASP_Setup module.
        energy_list = []
        energy_shift_dict = {"V": 1.682, "Cr": 2.013, "Mn": 1.681, "Fe": 2.733, "Co": 1.874, "Ni": 2.164, "O": 0.7023}

        if self.get_data_from_VASP_files == bool(True):
            pa = PoscarAnalyzer(poscar=self.poscar)
            oa = OszicarAnalyzer(oszicar=self.oszicar, poscar=self.poscar)
            energy = oa.get_energy()
            composition_dict = pa.get_composition_dict()
            print "The unshifted energy per cell is:", energy
            # Find which elements in material correspond to those that need shifting, if any, and apply energy shift
            for key in composition_dict.keys():
                if key in energy_shift_dict.keys():
                    energy -= energy_shift_dict[key]*composition_dict[key]
            print "The shifted energy per cell is:", energy
            energy_list.append(energy)

        elif self.get_data_from_VASP_files == bool(False):
            if len(self.composition_energy) > 0:
                for index, entry in enumerate(self.composition_energy):
                    energy = entry
                    print "The unshifted energy is:", energy

                    # Find which elements in material correspond to those that need shifting, if any, and apply energy shift
                    for key in self.composition_dict[index].keys():
                        if key in energy_shift_dict.keys():
                            energy -= energy_shift_dict[key]*self.composition_dict[index][key]
                    print "The shifted energy is:", energy
                    energy_list.append(energy)

        return energy_list

    def _create_chemical_system_from_MP(self):
        mp = MPRester(self.mapi_key)
        pa = PoscarAnalyzer(poscar=self.poscar)
        element_names = pa.get_element_names()

        # Remove redundant elements from element_names list
        elements_reduced = []
        for entry in element_names:
            if entry not in elements_reduced:
                elements_reduced.append(entry)
        # Remake elements list using only the reduced list of elements
        element_names = []
        for entry in elements_reduced:
                element_names.append(entry)

        if self.additional_elements_to_include is not None:
            if len(self.additional_elements_to_include) > 0:
                for entry in self.additional_elements_to_include:
                    if entry not in element_names:
                        print "%s is being added to the list of elements ..." % str(entry)
                        element_names.append(entry)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            print "Getting materials info from MP for the elements %s" % element_names, ", this may take a minute ..."
            entries_in_system = mp.get_entries_in_chemsys(element_names)

        return entries_in_system

    def _create_pdentry(self):
        pd_entry_list = []
        energy_list = self._calc_GGA_GGAU_shifted_energy()
        if self.get_data_from_VASP_files == bool(True):
            pa = PoscarAnalyzer(poscar=self.poscar)
            composition_dict = pa.get_composition_dict()
            comp = Composition(composition_dict)
            pd_entry = PDEntry(comp, energy_list[0])
            pd_entry_list.append(pd_entry)
        elif self.get_data_from_VASP_files == bool(False):
            for index, composition in enumerate(self.composition_dict):
                comp = Composition(composition)
                pd_entry = PDEntry(comp, energy_list[index])
                print "The PDEntry given for this system is:", pd_entry, " eV/cell"
                pd_entry_list.append(pd_entry)
        return pd_entry_list

    def _create_chemical_system_from_MP_and_remove_endmembers(self, species_to_remove):
        """
        Function which removes the entries from the MP chemical system so that the user can manually set the
        values of the chemical potentials. For example, if the user is looking at the La-Sr-Co-O system, but doesn't
        want to use the DFT value of the O chemical potential (which is unrealistically oxidizing), then the user can
        have the O endmember entries removed and insert their own O chemical potential value

        args:
            species_to_remove: list of elements to remove from entries list, e.g. ["O", "H"] will remove all O and H
                                end members so that the user can specify their own end member energies (for instance,
                                as a realistic chemical potential value). Viable elements to remove are:
                                [O, H, N, F, Cl, Br, I]
        """
        entries_in_system = self._create_chemical_system_from_MP()
        pd_entry = self._create_pdentry()

        # Remove the MP entries for elements specified by user. Need to loop multiple times over all entries in the
        # chemical system so that all desired entries are removed
        mp_list_O = ["mp-610917", "mp-607540", "mp-12957", "mp-973916", "mp-560602", "mp-611836"]
        mp_list_H = ["mp-632250", "mp-973783", "mp-1001570", "mp-24504", "mp-632172", "mp-632291", "mp-754417", "mp-850274", "mp-570752", "mp-23907", "mp-634659"]
        mp_list_N = ["mp-25", "mp-154", "mp-568584", "mp-999498", "mp-570747", "mp-12103", "mp-672234", "mp-672233", "mp-754514"]
        mp_list_F = ["mp-561367", "mp-21848", "mp-561203"]
        mp_list_Cl = ["mp-570778", "mp-22848"]
        mp_list_Br = ["mp-23154", "mp-998864", "mp-673171", "mp-998861"]
        mp_list_I = ["mp-23153", "mp-639751", "mp-601148", "mp-684663"]
        mp_list_together = [mp_list_O, mp_list_H, mp_list_N, mp_list_F, mp_list_Cl, mp_list_Br, mp_list_I]
        mp_dict = {"O": 0, "H": 1, "N": 2, "F": 3, "Cl": 4, "Br": 5, "I": 6}

        for species in species_to_remove:
            count = 0
            mp_list_to_use = mp_list_together[mp_dict[species]]
            while count < len(mp_list_to_use):
                for entry in entries_in_system:
                    if entry != pd_entry:
                        if entry.entry_id in mp_list_to_use:
                            print "removing an %s entry," % species, " entry ID number", entry.entry_id
                            entries_in_system.remove(entry)
                            count += 1

        # Remove material compositions specified by the user
        if self.material_ids_to_remove is not None:
            if len(self.material_ids_to_remove) > 0:
                for entry in self.material_ids_to_remove:
                    for mpentry in entries_in_system:
                        if mpentry.entry_id == entry:
                            print "removing material %s" % mpentry.entry_id
                            entries_in_system.remove(mpentry)

        return entries_in_system

    def _add_organic_molecules(self, add_roomtemp_gas_shift=False):
        element_names = PoscarAnalyzer(self.poscar).get_element_names()
        list_of_organic_compositions = []
        list_of_organic_energies = []

        # Can only use organic entries if C is an element in the POSCAR. Further, can only use N-containing organics if
        # N is an element in POSCAR, and only use halogen-containing organics if a halogen is in the POSCAR
        if "C" in element_names:
            organic_compositions_C = ['C2H2', 'C2H6', 'C3H6', 'C3H6', 'C4H2', 'C4H4', 'C4H6', 'C4H6', 'C4H6', 'C4H8', 'C4H8', 'C4H8', 'CH4'] #13 entries
            if add_roomtemp_gas_shift == False:
                organic_energies_C = [-22.9470, -40.4940, -45.4700, -48.4310, -34.1550, -44.3100, -56.6340, -56.5630, -57.0170, -65.3350, -65.1770, -65.3090, -24.0310] #13 entries
            if add_roomtemp_gas_shift == True:
                organic_energies_C = [-23.0720, -39.3950, -44.4480, -47.1940, -34.2200, -43.6910, -55.3830, -55.2210, -55.7180, -63.4380, -63.2880, -63.4300, -23.6790]
            list_of_organic_compositions.append(organic_compositions_C)
            list_of_organic_energies.append(organic_energies_C)

            if "N" in element_names:
                organic_compositions_C_N = ['C2H3N', 'C2H3N', 'C2H3N', 'C2H3N', 'C2H4N', 'C2H5N', 'C2H5N', 'C2H7N', 'C2HN', 'C3H3N', 'C3H3N', 'C3H3N', 'C3H3N', 'C3H5N', 'C3H5N', 'C3H5N', 'C3H5N', 'C3H5N', 'C3H5N', 'C3H6N', 'C3H7N', 'C3H7N', 'C3H7N', 'C3H7N', 'C3H7N', 'H3N'] #27 entries
                if add_roomtemp_gas_shift == False:
                    organic_energies_C_N = [-35.016, -26.716, -35.839, -33.447, -31.971, -44.195, -44.183, -51.924, -19.692, -43.049, -44.818, -43.876, -41.658, -52.354, -51.558, -51.710, -52.500, -51.872, -51.229, -45.470, -60.783, -60.139, -60.574, -60.286, -60.305, -19.520] #27 entries
                if add_roomtemp_gas_shift == True:
                    organic_energies_C_N = [-34.772, -26.481, -35.638, -33.261, -31.485, -43.324, -43.330, -50.435, -20.126, -42.707, -44.482, -43.555, -41.325, -51.376, -50.531, -50.791, -51.533, -50.842, -50.214, -44.425, -59.151, -58.716, -59.019, -58.634, -58.711, -19.465]
                list_of_organic_compositions.append(organic_compositions_C_N)
                list_of_organic_energies.append(organic_energies_C_N)

            if "I" in element_names:
                organic_compositions_C_I = ['C2H5I', 'C2H6I', 'C3H5I', 'C3H5I', 'C3H5I', 'C4H5I', 'C4H5I', 'C4H7I', 'C4H7I', 'C4H7I', 'C4H7I', 'C4H7I', 'C4H7I', 'C4H9I', 'CH3I'] #15 entries
                if add_roomtemp_gas_shift == False:
                    organic_energies_C_I = [-37.8650, -39.4240, -46.0820, -46.0580, -45.7520, -49.5710, -53.9720, -62.5650, -62.3980, -62.7410, -62.6740, -54.4520, -62.5350, -70.8140, -21.2850] #15 entries
                if add_roomtemp_gas_shift == True:
                    organic_energies_C_I = [-37.174482611, -38.8358178574, -45.3108410214, -45.2915944392, -44.9167679226, -48.8134780064, -53.1403520958, -60.790483545, -61.2510331184, -61.061126093, -61.1843491184, -53.0387394042, -61.0319775324, -68.6714105754, -21.2429180782]
                list_of_organic_compositions.append(organic_compositions_C_I)
                list_of_organic_energies.append(organic_energies_C_I)

            if "I" in element_names and "N" in element_names:
                organic_compositions_C_I_N = ['C2H6IN', 'C3H8IN', 'C3IN', 'CH3NH3I'] #4 entries
                if add_roomtemp_gas_shift == False:
                    organic_energies_C_I_N = [-49.6710, -65.8640, -33.2160, -41.1197] #4 entries
                if add_roomtemp_gas_shift == True:
                    organic_energies_C_I_N = [-48.572, -64.050, -33.838, -41.1197] #4 entries
                list_of_organic_compositions.append(organic_compositions_C_I_N)
                list_of_organic_energies.append(organic_energies_C_I_N)

            if "Br" in element_names:
                organic_compositions_C_Br = ['C2H5Br', 'C2H6Br', 'C3H5Br', 'C3H5Br', 'C3H5Br', 'C4H5Br', 'C4H5Br', 'C4H7Br', 'C4H7Br', 'C4H7Br', 'C4H7Br', 'C4H7Br', 'C4H7Br', 'C4H9Br', 'CH3Br'] #15 entries
                if add_roomtemp_gas_shift == False:
                    organic_energies_C_Br = [-38.409, -39.736, -46.662, -46.607, -46.297, -50.112, -54.468, -63.160, -62.964, -63.256, -63.185, -55.027, -63.071, -71.391, -21.790] #15 entries
                if add_roomtemp_gas_shift == True:
                    organic_energies_C_Br = [-37.6854620284, -39.0459850504, -45.859333928, -45.8093556546, -45.4304516416, -49.5406921368, -53.6034491474, -61.6264191706, -61.3244463448, -61.7313399492, -61.6602841606, -53.5816926434, -61.534848584, -69.2175999048, -21.7948053608]
                list_of_organic_compositions.append(organic_compositions_C_Br)
                list_of_organic_energies.append(organic_energies_C_Br)

            if "Br" in element_names and "N" in element_names:
                organic_compositions_C_Br_N = ['C2H6BrN', 'C3H8BrN', 'C3BrN'] #3 entries
                if add_roomtemp_gas_shift == False:
                    organic_energies_C_Br_N = [-50.189, -66.393, -33.630] #3 entries
                if add_roomtemp_gas_shift == True:
                    organic_energies_C_Br_N = [-49.061, -64.543, -34.221] #3 entries
                list_of_organic_compositions.append(organic_compositions_C_Br_N)
                list_of_organic_energies.append(organic_energies_C_Br_N)

            if "Cl" in element_names:
                organic_compositions_C_Cl = ['C2H5Cl', 'C2H6Cl', 'C3H5Cl', 'C3H5Cl', 'C3H5Cl', 'C4H5Cl', 'C4H5Cl', 'C4H7Cl', 'C4H7Cl', 'C4H7Cl', 'C4H7Cl', 'C4H7Cl', 'C4H7Cl', 'C4H9Cl', 'CH3Cl'] #15 entries
                if add_roomtemp_gas_shift == False:
                    organic_energies_C_Cl = [-39.038, -40.095, -47.352, -47.292, -46.950, -50.808, -55.096, -63.832, -63.609, -63.850, -63.775, -55.664, -63.700, -72.028, -22.409] #15 entries
                if add_roomtemp_gas_shift == True:
                    organic_energies_C_Cl = [-38.2659866224, -39.4707678498, -46.5060786144, -46.449482056, -46.0393081736, -50.2156424966, -54.1843084492, -62.2522185554, -61.9254377224, -62.2798507568, -62.2052597568, -54.173180451, -62.1158036744, -69.8104274202, -22.366939444]
                list_of_organic_compositions.append(organic_compositions_C_Cl)
                list_of_organic_energies.append(organic_energies_C_Cl)

            if "Cl" in element_names and "N" in element_names:
                organic_compositions_C_Cl_N = ['C2H6ClN', 'C3H8ClN', 'C3ClN'] #3 entries
                if add_roomtemp_gas_shift == False:
                    organic_energies_C_Cl_N = [-50.804, -67.022, -34.336] #3 entries
                if add_roomtemp_gas_shift == True:
                    organic_energies_C_Cl_N = [-49.631, -65.122, -34.882] #3 entries
                list_of_organic_compositions.append(organic_compositions_C_Cl_N)
                list_of_organic_energies.append(organic_energies_C_Cl_N)

            if "F" in element_names:
                organic_compositions_C_F = ['C2H5F', 'C2H6F', 'C3H5F', 'C3H5F', 'C3H5F', 'C4H5F', 'C4H5F', 'C4H7F', 'C4H7F', 'C4H7F', 'C4H7F', 'C4H7F', 'C4H7F', 'C4H9F', 'CH3F'] #15 entries
                if add_roomtemp_gas_shift == False:
                    organic_energies_C_F = [-40.681, -42.113, -49.104, -47.948, -48.577, -50.216, -56.036, -65.602, -65.309, -65.465, -65.388, -57.389, -65.330, -73.776, -21.790] #15 entries
                if add_roomtemp_gas_shift == True:
                    organic_energies_C_F = [-39.8477529154, -41.2884478318, -48.1936499218, -47.0411676152, -47.603645279, -49.447237007, -55.0605830012, -63.9604352784, -63.5608601272, -63.8323205772, -63.7549077886, -55.8429766748, -63.6839136088, -71.5002315776, -21.683864483]
                list_of_organic_compositions.append(organic_compositions_C_F)
                list_of_organic_energies.append(organic_energies_C_F)

            if "F" in element_names and "N" in element_names:
                organic_compositions_C_F_N = ['C2H6FN', 'C3H8FN', 'C3FN'] #3 entries
                if add_roomtemp_gas_shift == False:
                    organic_energies_C_F_N = [-52.395, -68.629, -35.443] #3 entries
                if add_roomtemp_gas_shift == True:
                    organic_energies_C_F_N = [-51.161, -66.663, -35.925] #3 entries
                list_of_organic_compositions.append(organic_compositions_C_F_N)
                list_of_organic_energies.append(organic_energies_C_F_N)

        organic_compositions = []
        organic_energies = []
        for index in range(len(list_of_organic_compositions)):
            for entry1, entry2 in zip(list_of_organic_compositions[index], list_of_organic_energies[index]):
                organic_compositions.append(entry1)
                organic_energies.append(entry2)

        print "Adding %i organic molecules" % len(organic_compositions)
        print "The list of organic molecules is:"
        print organic_compositions
        print "The list of organic molecule energies is:"
        print organic_energies

        return organic_compositions, organic_energies

class ChemicalPotentialAnalyzer():
    """Class used to calculate the chemical potential of select gaseous species under user-specified environmental conditions
    The elements that can have their chemical potential calculated are: O, H, N, F, Cl, Br, I
    args:
        temperature: (int) temperature in K, must be > 0
        pressure: (int) pressure in atm, must be > 0
        functional: (str) DFT XC functional, specify as "PBE" or "PW91". Note that getting mu_H from mu_H2O only works
                    for PBE.
        relative_humidity: (float) number between 0 and 1, specifies the % relative humidity of H2O. Only used if you are
            acquiring the H chemical potential via equilibrium with H2O.
        energy_shift: (bool) whether to add the energy correction to O energy. Note that for Stability calculations,
            use shift of "False" because MP convention is to apply shift to compounds, not gaseous O.
        relative_humidity : (float) value of relative humidity, which determines partial pressure of H2O
        pressure_O2_forH2O: (float) value for partial pressure of O2, only for use in getting chem pot of H from H2O. Defaults to
                    0.2 atm, which is the partial pressure of O2 in air.
    instance methods:
        get_O_chem_pot : (float) calculates the O gas chemical potential, enthalpy and entropy
        get_H_chem_pot : (float) calculates the H gas chemical potential, enthalpy and entropy
        get_H_chem_pot_fromH2O : (float) calculates the H chemical potential assuming that H2 is in equilibrium with H2O
        get_N_chem_pot : (float) calculates the N gas chemical potential, enthalpy and entropy
        get_F_chem_pot : (float) calculates the F gas chemical potential, enthalpy and entropy
        get_Cl_chem_pot : (float) calculates the Cl gas chemical potential, enthalpy and entropy
        get_Br_chem_pot : (float) calculates the Br gas chemical potential, enthalpy and entropy
        get_I_chem_pot : (float) calculates the I gas chemical potential, enthalpy and entropy
    """
    def __init__(self, temperature, pressure, functional, energy_shift=bool(False), relative_humidity=1, pressure_O2_forH2O=0.2):
        self.temperature = temperature
        self.pressure = pressure
        self.functional = functional
        self.energy_shift = energy_shift
        self.relative_humidity = relative_humidity
        self.pressure_O2_forH2O = pressure_O2_forH2O

    def get_O_chem_pot(self):
        self._check_input()
        k = 8.61733*10**-5 #Boltzmann's constant in eV/K
        P0 = 1 #standard pressure in atm
        O_shift_PBE = 2*0.7023
        O_shift_PW91 = 2*0.165
        conv = 1.60217662*10**-19 #J per eV
        Nav = 6.0221409*10**23 #Avogadro's number
        EO2_PBE = -9.871 #Materials Project value of O2 energy, PBE
        EO2_PW91 = -9.09 #Yueh-Lin Lee PRB 2009 value for O2 energy, PW91

        t = float(float(self.temperature)/1000)
        # These constants good for T = 100-700
        A1_O2 = 31.32234
        B1_O2 = -20.23531
        C1_O2 = 57.86644
        D1_O2 = -36.50624
        E1_O2 = -0.007374
        F1_O2 = -8.903471
        G1_O2 = 246.7945
        # These constants good for T = 700-2000
        A2_O2 = 30.03235
        B2_O2 = 8.772972
        C2_O2 = -3.988133
        D2_O2 = 0.788313
        E2_O2 = -0.741599
        F2_O2 = -11.32468
        G2_O2 = 236.1663

        if self.temperature <= 700:
            S0 = A1_O2*math.log(t) + B1_O2*t + (C1_O2*t**2)/2 + (D1_O2*t**3)/3 - E1_O2/(2*t**2) + G1_O2  #in J/mol-K
            HT_H0 = A1_O2*t + (B1_O2*t**2)/2 + (C1_O2*t**3)/3 + (D1_O2*t**4)/4 - E1_O2/(t) + F1_O2  #in kJ/mol
            S0 = S0 /(conv*Nav) #in eV/atom-K
            HT_H0 = (1000*HT_H0) / (conv*Nav) #in eV/atom

            if self.functional == "PBE" or self.functional == "pbe":
                if self.energy_shift == False:
                    mu = (EO2_PBE + HT_H0 - self.temperature*S0 + k*self.temperature*math.log(self.pressure/P0))/2
                if self.energy_shift == True:
                    mu = (EO2_PBE + O_shift_PBE + HT_H0 - self.temperature*S0 + k*self.temperature*math.log(self.pressure/P0))/2

            if self.functional == "PW91" or self.functional == "pw91":
                if self.energy_shift == False:
                    mu = (EO2_PW91 + HT_H0 - self.temperature*S0 + k*self.temperature*math.log(self.pressure/P0))/2
                if self.energy_shift == True:
                    mu = (EO2_PW91 + O_shift_PW91 + HT_H0 - self.temperature*S0 + k*self.temperature*math.log(self.pressure/P0))/2

        if self.temperature > 700:
            S0 = A2_O2*math.log(t) + B2_O2*t + (C2_O2*t**2)/2 + (D2_O2*t**3)/3 - E2_O2/(2*t**2) + G2_O2  #in J/mol-K
            HT_H0 = A2_O2*t + (B2_O2*t**2)/2 + (C2_O2*t**3)/3 + (D2_O2*t**4)/4 - E2_O2/(t) + F2_O2  #in kJ/mol
            S0 = S0 /(conv*Nav) #in eV/atom-K
            HT_H0 = (1000*HT_H0) / (conv*Nav) #in eV/atom

            if self.functional == "PBE" or self.functional == "pbe":
                if self.energy_shift == False:
                    mu = (EO2_PBE + HT_H0 - self.temperature*S0 + k*self.temperature*math.log(self.pressure/P0))/2
                if self.energy_shift == True:
                    mu = (EO2_PBE + O_shift_PBE + HT_H0 - self.temperature*S0 + k*self.temperature*math.log(self.pressure/P0))/2

            if self.functional == "PW91" or self.functional == "pw91":
                if self.energy_shift == False:
                    mu = (EO2_PW91 + HT_H0 - self.temperature*S0 + k*self.temperature*math.log(self.pressure/P0))/2
                if self.energy_shift == True:
                    mu = (EO2_PW91 + O_shift_PW91 + HT_H0 - self.temperature*S0 + k*self.temperature*math.log(self.pressure/P0))/2
        return mu, HT_H0, S0

    def get_H_chem_pot(self):
        self._check_input()
        k = 8.61733*10**-5 #Boltzmann's constant in eV/K
        P0 = 1 #standard pressure in atm
        H2_shift_from0K = -0.087766 #Shift, in eV/H2, of enthalpy of H2 from T=298K to T=0K. Subtract this so that DFT energy more comparable to T=298K standard state
        ZPE_H2 = 0.2699 #eV per H2, zero-point energy. This value obtained from wavenumber of 2179.3 1/cm, from Experimental Vibrational Zero-Point Energies: Diatomic Molecules Karl K. Irikuraa Physical and Chemical Properties Division, National Institute of Standards and Technology, Gaithersburg, Maryland 20899-83
        conv = 1.60217662*10**-19 #J per eV
        Nav = 6.0221409*10**23 #Avogadro's number
        EH2_PBE = -6.717 #DFT energy of H2, PBE in eV/H2
        EH2_PW91 = -6.746 #DFT energy of H2, PW91 in eV/H2

        t = float(float(self.temperature)/1000)
        # These constants good for T = 298-1000
        A1_H2 = 33.066178
        B1_H2 = -11.363417
        C1_H2 = 11.432816
        D1_H2 = -2.772874
        E1_H2 = -0.158558
        F1_H2 = -9.980797
        G1_H2 = 172.707974

        # These constants good for T = 1000-2500
        A2_H2 = 18.563083
        B2_H2 = 12.257357
        C2_H2 = -2.859786
        D2_H2 = 0.268238
        E2_H2 = 1.977990
        F2_H2 = -1.147438
        G2_H2 = 156.288133

        if self.temperature <= 1000:
            S0 = A1_H2*math.log(t) + B1_H2*t + (C1_H2*t**2)/2 + (D1_H2*t**3)/3 - E1_H2/(2*t**2) + G1_H2  #in J/mol-K
            HT_H0 = A1_H2*t + (B1_H2*t**2)/2 + (C1_H2*t**3)/3 + (D1_H2*t**4)/4 - E1_H2/(t) + F1_H2  #in kJ/mol
            S0 = S0 /(conv*Nav) #in eV/atom-K
            HT_H0 = (1000*HT_H0) / (conv*Nav) #in eV/atom

            if self.functional == "PBE" or self.functional == "pbe":
                if self.energy_shift == False:
                    mu = (EH2_PBE - H2_shift_from0K + HT_H0 - self.temperature*S0 + k*self.temperature*math.log(self.pressure/P0))/2
                if self.energy_shift == True:
                    mu = (EH2_PBE - H2_shift_from0K + ZPE_H2 + HT_H0 - self.temperature*S0 + k*self.temperature*math.log(self.pressure/P0))/2

            if self.functional == "PW91" or self.functional == "pw91":
                if self.energy_shift == False:
                    mu = (EH2_PW91 - H2_shift_from0K + HT_H0 - self.temperature*S0 + k*self.temperature*math.log(self.pressure/P0))/2
                if self.energy_shift == True:
                    mu = (EH2_PW91 - H2_shift_from0K + ZPE_H2 + HT_H0 - self.temperature*S0 + k*self.temperature*math.log(self.pressure/P0))/2

        if self.temperature > 1000:
            S0 = A2_H2*math.log(t) + B2_H2*t + (C2_H2*t**2)/2 + (D2_H2*t**3)/3 - E2_H2/(2*t**2) + G2_H2  #in J/mol-K
            HT_H0 = A2_H2*t + (B2_H2*t**2)/2 + (C2_H2*t**3)/3 + (D2_H2*t**4)/4 - E2_H2/(t) + F2_H2  #in kJ/mol
            S0 = S0 /(conv*Nav) #in eV/atom-K
            HT_H0 = (1000*HT_H0) / (conv*Nav) #in eV/atom

            if self.functional == "PBE" or self.functional == "pbe":
                if self.energy_shift == False:
                    mu = (EH2_PBE - H2_shift_from0K + HT_H0 - self.temperature*S0 + k*self.temperature*math.log(self.pressure/P0))/2
                if self.energy_shift == True:
                    mu = (EH2_PBE - H2_shift_from0K + ZPE_H2 + HT_H0 - self.temperature*S0 + k*self.temperature*math.log(self.pressure/P0))/2

            if self.functional == "PW91" or self.functional == "pw91":
                if self.energy_shift == False:
                    mu = (EH2_PW91 - H2_shift_from0K + HT_H0 - self.temperature*S0 + k*self.temperature*math.log(self.pressure/P0))/2
                if self.energy_shift == True:
                    mu = (EH2_PW91 - H2_shift_from0K + ZPE_H2 + HT_H0 - self.temperature*S0 + k*self.temperature*math.log(self.pressure/P0))/2

        return mu, HT_H0, S0

    def get_H_chem_pot_fromH2O(self):
        # A note of caution when using this function: The H2O energies used here are ONLY for PBE. A similar analysis
        # with PW91 has not been performed.
        if self.functional != "PBE":
            raise TypeError('This function is only applicable to DFT-PBE energies, exiting...')
        # A conceptual note about why the chemical potential of H2O at 298 is used: The O and H chemical potentials are
        # changed, so H2O DFT energy remains fixed (only end members change). This choice of chemical potentials reproduces
        # the correct experimental formation energy of H2O at the relevant conditions of interest to the user. This is
        # how the phase diagrams in MP are calculated: only end member energy changes, and this changes formation energies
        # of all compounds, however the effective calculated energy of each compound is unchanged
        k = 8.61733*10**-5 #Boltzmann's constant in eV/K
        conv = 1.60217662*10**-19 #J per eV
        Nav = 6.0221409*10**23 #Avogadro's number k = 8.61733*10**-5 #Boltzmann's const
        P0 = 1 #standard pressure in atm
        H_melt = 0.0623 #eV/H2O, at T=273K
        H_vaporization_298 = 0.4568 #eV/H2O, at T=298K
        H_vaporization_373 = 0.4222 #eV/H2O, at T=373K
        Hform0_liqH2O = -2.9665 #eV/H2O standard enthalpy of formation for liquid at 298K
        Hform0_gasH2O = -2.5098 #eV/H2O standard enthalpy of formation for gas at 298K

        # Getting partial pressure of H2O from specified relative humidity using Antoine equation

        A = 5.40221 ; B = 1838.675 ; C = -31.737 ; T = 298
        # Equilibrium vapor pressure of H2O at 298 K:
        P_H2O = 10**(A-(B/(C+T)))  #atm
        # Adjust pressure for relative humidity
        P_H2O = P_H2O*self.relative_humidity #atm
        mu_O, H_O, S_O = self.get_O_chem_pot()
        mu_H_gas, H_H, S_H = self.get_H_chem_pot()

        mu_O0, H0_O, S0_O = ChemicalPotentialAnalyzer(temperature=298, pressure=1, functional=self.functional,
                                                      relative_humidity=1, energy_shift=self.energy_shift).get_O_chem_pot()
        mu_H0_gas, H0_H, S0_H = ChemicalPotentialAnalyzer(temperature=298, pressure=1, functional=self.functional,
                                                          relative_humidity=1, energy_shift=self.energy_shift).get_H_chem_pot()

        t = float(float(self.temperature)/1000)
        # These constants good for T = 298-500 (H2O in liquid phase)
        A1_H2O = -203.606
        B1_H2O = 1523.29
        C1_H2O = -3196.413
        D1_H2O = 2474.455
        E1_H2O = 3.855326
        F1_H2O = -256.5478
        G1_H2O = -488.7163
        H1_H2O = -285.8304

        # These constants good for T = 500-1700 (H2O in gas phase). USE THIS FOR GAS PHASE ENERGIES
        A2_H2O = 30.092
        B2_H2O = 6.832514
        C2_H2O = 6.793435
        D2_H2O = -2.53448
        E2_H2O = 0.082139
        F2_H2O = -250.881
        G2_H2O = 223.3967
        H2_H2O = -241.8264

        S_H2O = A2_H2O*math.log(t) + B2_H2O*(t) + (C2_H2O*(t)**2)/2 + (D2_H2O*(t)**3)/3 - E2_H2O/(2*(t)**2) + G2_H2O
        H_H2O = A2_H2O*(t) + (B2_H2O*(t)**2)/2 + (C2_H2O*(t)**3)/3 + (D2_H2O*(t)**4)/4 - E2_H2O/(t) + F2_H2O - H2_H2O

        t_const = 0.298 #reduced temperature 298/1000
        S0_H2O = A2_H2O*math.log(t_const) + B2_H2O*t_const + (C2_H2O*t_const**2)/2 + (D2_H2O*t_const**3)/3 - E2_H2O/(2*t_const**2) + G2_H2O #in J/mol-K
        H0_H2O = A2_H2O*t_const + (B2_H2O*t_const**2)/2 + (C2_H2O*t_const**3)/3 + (D2_H2O*t_const**4)/4 - E2_H2O/(t_const) + F2_H2O - H2_H2O #in kJ/mol

        S_H2O = S_H2O /(conv*Nav)  #in eV/atom-K
        H_H2O = (1000*H_H2O) / (conv*Nav) #in eV/atom
        S0_H2O = S0_H2O /(conv*Nav)  #in eV/atom-K
        H0_H2O = (1000*H0_H2O) / (conv*Nav) #in eV/atom

        Gform_H2O = (Hform0_liqH2O + H_vaporization_298 + H_H2O - H_H - 0.5*H_O - (H0_H2O - H0_H - 0.5*H0_O)) \
                    - self.temperature*(S_H2O - S_H - 0.5*S_O) -(S0_H2O- S0_H - 0.5*S0_O)*(self.temperature-298) \
                    + k*self.temperature*math.log(P_H2O/P0)

        # Need to incorporate O correction shift. If applying shift to O, then don't apply to H2O. If not on O, needs
        # to be on H2O.
        if self.energy_shift == bool(True):
            mu_H2O = -14.885 # MP calculated value, eV/H2O, no O shift because shift is applied to O gas, not compound
        if self.energy_shift == bool(False):
            mu_H2O = -14.885 - 0.7023 # Apply the O shift to H2O value, this is needed if doing stability calculations
            if self.functional == "PW91" or self.functional == "pw91":
                raise TypeError("This analysis with H2O is only valid for PBE functionals, any results with a PW91 O shift may be incorrect")

        mu_H = 0.5*(mu_H2O - mu_O - Gform_H2O)
        return mu_H

    def get_N_chem_pot(self):
        self._check_input()
        #Check to make sure PBE functional is specifed, else quit
        if self.functional != "PBE":
            raise TypeError('This function is only applicable to DFT-PBE energies, exiting...')
        k = 8.61733*10**-5 #Boltzmann's constant in eV/K
        P0 = 1 #standard pressure in atm
        conv = 1.60217662*10**-19 #J per eV
        Nav = 6.0221409*10**23 #Avogadro's number
        EN2_PBE = -16.653 #Materials Project value of N2 energy, PBE
        N2_shift_from0K = -0.0900465 #Shift, in eV/N2, of enthalpy of N2 from T=298K to T=0K. Subtract this so that DFT energy more comparable to T=298K standard state
        N_shift_PBE = 2*0.37425 #Materials Project gas shift for N2, eV/N2

        t = float(float(self.temperature)/1000)
        # These constants good for T = 100-500
        A1 = 28.98641
        B1 = 1.853978
        C1 = -9.647458
        D1 = 16.63537
        E1 = 0.000117
        F1 = -8.671914
        G1 = 226.4168
        # These constants good for T = 500-2000
        A2 = 19.50583
        B2 = 19.88705
        C2 = -8.598535
        D2 = 1.369784
        E2 = 0.527601
        F2 = -4.935202
        G2 = 212.39

        if self.temperature <= 500:
            S0 = A1*math.log(t) + B1*t + (C1*t**2)/2 + (D1*t**3)/3 - E1/(2*t**2) + G1  #in J/mol-K
            HT_H0 = A1*t + (B1*t**2)/2 + (C1*t**3)/3 + (D1*t**4)/4 - E1/(t) + F1  #in kJ/mol
            S0 = S0 /(conv*Nav) #in eV/atom-K
            HT_H0 = (1000*HT_H0) / (conv*Nav) #in eV/atom

            if self.functional == "PBE" or self.functional == "pbe":
                if self.energy_shift == False:
                    mu = (EN2_PBE - N2_shift_from0K + HT_H0 - self.temperature*S0 + k*self.temperature*math.log(self.pressure/P0))/2
                if self.energy_shift == True:
                    print "WARNING: You are calculating N2 energy with the materials project energy shift, this has not been well-tested!!!"
                    mu = (EN2_PBE + N_shift_PBE + HT_H0 - self.temperature*S0 + k*self.temperature*math.log(self.pressure/P0))/2

        if self.temperature > 500:
            S0 = A2*math.log(t) + B2*t + (C2*t**2)/2 + (D2*t**3)/3 - E2/(2*t**2) + G2  #in J/mol-K
            HT_H0 = A2*t + (B2*t**2)/2 + (C2*t**3)/3 + (D2*t**4)/4 - E2/(t) + F2  #in kJ/mol
            S0 = S0 /(conv*Nav) #in eV/atom-K
            HT_H0 = (1000*HT_H0) / (conv*Nav) #in eV/atom

            if self.functional == "PBE" or self.functional == "pbe":
                if self.energy_shift == False:
                    mu = (EN2_PBE - N2_shift_from0K + HT_H0 - self.temperature*S0 + k*self.temperature*math.log(self.pressure/P0))/2
                if self.energy_shift == True:
                    print "WARNING: You are calculating N2 energy with the materials project energy shift, this has not been well-tested!!!"
                    mu = (EN2_PBE + N_shift_PBE + HT_H0 - self.temperature*S0 + k*self.temperature*math.log(self.pressure/P0))/2

        return mu, HT_H0, S0

    def get_F_chem_pot(self):
        self._check_input()
        #Check to make sure PBE functional is specifed, else quit
        if self.functional != "PBE":
            raise TypeError('This function is only applicable to DFT-PBE energies, exiting...')
        k = 8.61733*10**-5 #Boltzmann's constant in eV/K
        P0 = 1 #standard pressure in atm
        conv = 1.60217662*10**-19 #J per eV
        Nav = 6.0221409*10**23 #Avogadro's number
        EF2_PBE = -3.7453 #Materials Project value of F2 energy, PBE
        F_shift_PBE = 2*0.4521 #Materials Project gas shift for F2, eV/F2
        F2_shift_from0K = -0.096947 #Shift, in eV/F2, of enthalpy of F2 from T=298K to T=0K. Subtract this so that DFT energy more comparable to T=298K standard state

        # Data for F starts here!

        t = float(float(self.temperature)/1000)
        # These constants good for T = 298-6000
        A1 = 31.4451
        B1 = 8.413831
        C1 = -2.77885
        D1 = 0.218194
        E1 = -0.211175
        F1 = -10.4326
        G1 = 237.277

        if self.temperature <= 6000:
            S0 = A1*math.log(t) + B1*t + (C1*t**2)/2 + (D1*t**3)/3 - E1/(2*t**2) + G1  #in J/mol-K
            HT_H0 = A1*t + (B1*t**2)/2 + (C1*t**3)/3 + (D1*t**4)/4 - E1/(t) + F1  #in kJ/mol
            S0 = S0 /(conv*Nav) #in eV/atom-K
            HT_H0 = (1000*HT_H0) / (conv*Nav) #in eV/atom

            if self.functional == "PBE" or self.functional == "pbe":
                if self.energy_shift == False:
                    mu = (EF2_PBE - F2_shift_from0K + HT_H0 - self.temperature*S0 + k*self.temperature*math.log(self.pressure/P0))/2
                if self.energy_shift == True:
                    print "WARNING: You are calculating F2 energy with the materials project energy shift, this has not been well-tested!!!"
                    mu = (EF2_PBE + F_shift_PBE + HT_H0 - self.temperature*S0 + k*self.temperature*math.log(self.pressure/P0))/2

        return mu, HT_H0, S0

    def get_Cl_chem_pot(self):
        self._check_input()
        #Check to make sure PBE functional is specifed, else quit
        if self.functional != "PBE":
            raise TypeError('This function is only applicable to DFT-PBE energies, exiting...')
        k = 8.61733*10**-5 #Boltzmann's constant in eV/K
        P0 = 1 #standard pressure in atm
        conv = 1.60217662*10**-19 #J per eV
        Nav = 6.0221409*10**23 #Avogadro's number
        ECl2_PBE = -3.5997 #Materials Project value of Cl2 energy, PBE
        Cl_shift_PBE = 2*0.5460 #Materials Project gas shift for Cl2, eV/Cl2
        Cl2_shift_from0K = -0.105023362 #Shift, in eV/Cl2, of enthalpy of Cl2 from T=298K to T=0K. Subtract this so that DFT energy more comparable to T=298K standard state

        t = float(float(self.temperature)/1000)
        # These constants good for T = 298-1000
        A1 = 33.0506
        B1 = 12.2294
        C1 = -12.0651
        D1 = 4.38533
        E1 = -0.159494
        F1 = -10.8348
        G1 = 259.029
        # These constants good for T > 1000
        A2 = 42.6773
        B2 = -5.00957
        C2 = 1.904621
        D2 = -0.165641
        E2 = -2.09848
        F2 = -17.2898
        G2 = 269.840

        if self.temperature <= 500:
            S0 = A1*math.log(t) + B1*t + (C1*t**2)/2 + (D1*t**3)/3 - E1/(2*t**2) + G1  #in J/mol-K
            HT_H0 = A1*t + (B1*t**2)/2 + (C1*t**3)/3 + (D1*t**4)/4 - E1/(t) + F1  #in kJ/mol
            S0 = S0 /(conv*Nav) #in eV/atom-K
            HT_H0 = (1000*HT_H0) / (conv*Nav) #in eV/atom

            if self.functional == "PBE" or self.functional == "pbe":
                if self.energy_shift == False:
                    mu = (ECl2_PBE - Cl2_shift_from0K + HT_H0 - self.temperature*S0 + k*self.temperature*math.log(self.pressure/P0))/2
                if self.energy_shift == True:
                    print "WARNING: You are calculating Cl2 energy with the materials project energy shift, this has not been well-tested!!!"
                    mu = (ECl2_PBE + Cl_shift_PBE + HT_H0 - self.temperature*S0 + k*self.temperature*math.log(self.pressure/P0))/2

        if self.temperature > 500:
            S0 = A2*math.log(t) + B2*t + (C2*t**2)/2 + (D2*t**3)/3 - E2/(2*t**2) + G2  #in J/mol-K
            HT_H0 = A2*t + (B2*t**2)/2 + (C2*t**3)/3 + (D2*t**4)/4 - E2/(t) + F2  #in kJ/mol
            S0 = S0 /(conv*Nav) #in eV/atom-K
            HT_H0 = (1000*HT_H0) / (conv*Nav) #in eV/atom

            if self.functional == "PBE" or self.functional == "pbe":
                if self.energy_shift == False:
                    mu = (ECl2_PBE - Cl2_shift_from0K + HT_H0 - self.temperature*S0 + k*self.temperature*math.log(self.pressure/P0))/2
                if self.energy_shift == True:
                    print "WARNING: You are calculating Cl2 energy with the materials project energy shift, this has not been well-tested!!!"
                    mu = (ECl2_PBE + Cl_shift_PBE + HT_H0 - self.temperature*S0 + k*self.temperature*math.log(self.pressure/P0))/2

        return mu, HT_H0, S0

    def get_Br_chem_pot(self):
        self._check_input()
        #Check to make sure PBE functional is specifed, else quit
        if self.functional != "PBE":
            raise TypeError('This function is only applicable to DFT-PBE energies, exiting...')
        k = 8.61733*10**-5 #Boltzmann's constant in eV/K
        P0 = 1 #standard pressure in atm
        conv = 1.60217662*10**-19 #J per eV
        Nav = 6.0221409*10**23 #Avogadro's number
        EBr2_PBE = -3.2549 #Materials Project value of Br2 energy, PBE
        Br_shift_PBE = 0 #Materials Project gas shift for Br2, eV/Br2
        Br2_shift_from0K = -0.112480826 #Shift, in eV/Br2, of enthalpy of Br2 from T=298K to T=0K. Subtract this so that DFT energy
        # more comparable to T=298K standard state. No liquid thermochem data for Br2 on NIST website, assume shift is for gaseous Br2.
        H_vap_Br2 = 0.324850548 #H of vaporization in eV/Br2. Br2 is liquid under standard conditions, and vaporizes at 332 K. DFT calculation
        # is for sol/liq Br2, so add this term to the free energy.

        t = float(float(self.temperature)/1000)
        # These constants good for T = 332-3400, Br2 in gaseous state. USE THIS FOR GAS PHASE ENERGIES
        A1 = 38.52723
        B1 = -1.976835
        C1 = 1.526107
        D1 = -0.198398
        E1 = -0.185815
        F1 = 18.8762
        G1 = 291.4863
        H1 = 30.91001

        if self.temperature <= 3400:
            S0 = A1*math.log(t) + B1*t + (C1*t**2)/2 + (D1*t**3)/3 - E1/(2*t**2) + G1  #in J/mol-K
            HT_H0 = A1*t + (B1*t**2)/2 + (C1*t**3)/3 + (D1*t**4)/4 - E1/(t) + F1 - H1  #in kJ/mol
            S0 = S0 /(conv*Nav) #in eV/atom-K
            HT_H0 = (1000*HT_H0) / (conv*Nav) #in eV/atom

            if self.functional == "PBE" or self.functional == "pbe":
                if self.energy_shift == False:
                    mu = (EBr2_PBE - Br2_shift_from0K + H_vap_Br2 + HT_H0 - self.temperature*S0 + k*self.temperature*math.log(self.pressure/P0))/2
                if self.energy_shift == True:
                    print "NOTE: There is currently no energy shift for Br2 implemented in Materials Project. You will get the same answer as if no shift is applied!"
                    mu = (EBr2_PBE - Br2_shift_from0K + H_vap_Br2 + Br_shift_PBE + HT_H0 - self.temperature*S0 + k*self.temperature*math.log(self.pressure/P0))/2

        return mu, HT_H0, S0

    def get_I_chem_pot(self):
        self._check_input()
        #Check to make sure PBE functional is specifed, else quit
        if self.functional != "PBE":
            raise TypeError('This function is only applicable to DFT-PBE energies, exiting...')
        k = 8.61733*10**-5 #Boltzmann's constant in eV/K
        P0 = 1 #standard pressure in atm
        conv = 1.60217662*10**-19 #J per eV
        Nav = 6.0221409*10**23 #Avogadro's number
        EI2_PBE = -3.0378 #Materials Project value of I2 energy, PBE
        I_shift_PBE = 0 #Materials Project gas shift for I2, eV/I2
        I2_shift_from0K = -0.115760872 #Shift, in eV/I2, of enthalpy of I2 from T=298K to T=0K. Subtract this so that DFT energy
        # more comparable to T=298K standard state. No solid thermochem data for I2 on NIST website, assume shift is for gaseous I2.
        H_fusion_I2 = 0.161076054 #H of melting in eV/I2. I2 is solid under standard conditions, and melts at 387 K. We want energy of I2 gas at
        # all temperatures, so add this term into free energy
        H_vap_I2 = 0.431750249 #H of vaporization in eV/I2. I2 is solid under standard conditions, and vaporizes at 458 K. We want energy of I2 gas at
        # all temperatures, so add this term into free energy

        t = float(float(self.temperature)/1000)
        # These constants good for T = 298-387 (I2 is solid)
        A1 = -195.7635
        B1 = 918.8984
        C1 = -1079.242
        D1 = 535.3219
        E1 = 5.156403
        F1 = 43.29938
        G1 = -322.478
        # These constants good for T > 387 and T < 458 (I2 is liquid)
        A2 = 80.66919
        B2 = 6.85*10**-8
        C2 = -8.72*10**-10
        D2 = 3.72313*10**-8
        E2 = 4.73583*10**-10
        F2 = -10.52782
        G2 = 247.9798
        H2 = 13.52302
        # These constants good for T > 458 (I2 is gas) USE THIS FOR GAS PHASE ENERGIES
        A3 = 37.79763
        B3 = 0.225453
        C3 = -0.912556
        D3 = 1.034913
        E3 = -0.083826
        F3 = 50.86865
        G3 = 305.9199
        H3 = 62.4211

        S0 = A3*math.log(t) + B3*t + (C3*t**2)/2 + (D3*t**3)/3 - E3/(2*t**2) + G3  #in J/mol-K
        HT_H0 = A3*t + (B3*t**2)/2 + (C3*t**3)/3 + (D3*t**4)/4 - E3/(t) + F3 - H3  #in kJ/mol
        S0 = S0 /(conv*Nav) #in eV/atom-K
        HT_H0 = (1000*HT_H0) / (conv*Nav) #in eV/atom

        if self.functional == "PBE" or self.functional == "pbe":
            if self.energy_shift == False:
                mu = (EI2_PBE - I2_shift_from0K + H_fusion_I2 + H_vap_I2 + HT_H0 - self.temperature*S0 + k*self.temperature*math.log(self.pressure/P0))/2
            if self.energy_shift == True:
                print "NOTE: There is currently no energy shift for I2 implemented in Materials Project. You will get the same answer as if no shift is applied!"
                mu = (EI2_PBE - I2_shift_from0K + H_fusion_I2 + H_vap_I2 + I_shift_PBE + HT_H0 - self.temperature*S0 + k*self.temperature*math.log(self.pressure/P0))/2

        return mu, HT_H0, S0

    def _check_input(self):
        func_list = ["pw91", "PW91", "PBE", "pbe"]
        if self.temperature <= 0:
            print "ERROR: You have selected on a nonsensical temperature. Exiting..."
            exit()
        if self.pressure <= 0:
            print "ERROR: You have selected on a nonsensical pressure. Exiting..."
            exit()
        if self.functional not in func_list:
            print "ERROR: You must specify one of either PW91 or PBE functionals. Exiting..."
            exit()
        # if self.energy_shift == "True":
        #    print "Applying the shift corrections"
        # if self.energy_shift == "False":
        #    print "No shift corrections made"
        return None
