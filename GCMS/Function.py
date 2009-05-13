"""
Provides generic IO functions for GC-MS data objects
"""

 #############################################################################
 #                                                                           #
 #    PyMS software for processing of metabolomic mass-spectrometry data     #
 #    Copyright (C) 2005-8 Vladimir Likic                                    #
 #                                                                           #
 #    This program is free software; you can redistribute it and/or modify   #
 #    it under the terms of the GNU General Public License version 2 as      #
 #    published by the Free Software Foundation.                             #
 #                                                                           #
 #    This program is distributed in the hope that it will be useful,        #
 #    but WITHOUT ANY WARRANTY; without even the implied warranty of         #
 #    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          #
 #    GNU General Public License for more details.                           #
 #                                                                           #
 #    You should have received a copy of the GNU General Public License      #
 #    along with this program; if not, write to the Free Software            #
 #    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.              #
 #                                                                           #
 #############################################################################

import math

from pyms.Utils.Error import error
from pyms.Utils.Utils import is_number, is_str, is_array, is_list, is_int
from pyms.GCMS.Class import GCMS_data, IntensityMatrix

def build_intensity_matrix(data, bin_size = 1):

    """
    @summary: Sets the full intensity matrix with flexible bins

    @param data: Raw GCMS data
    @type bin_size: pyms.GCMS.GCMS_data

    @param bin_size: bin size
    @type bin_size: IntType or FloatType

    @return: Binned IntensityMatrix object
    @rtype: pyms.GCMS.IntensityMatrix

    @author: Qiao Wang
    @author: Andrew Isaac
    @author: Vladimir Likic
    """

    if not isinstance(data, GCMS_data):
        error("data must be an GCMS_data object")
    if bin_size <= 0:
        error("The bin size must be larger than zero.")

    min_mass = data.get_min_mass()
    max_mass = data.get_max_mass()
    return __fill_bins(data, min_mass, max_mass, bin_size)

def build_intensity_matrix_i(data):

    """
    @summary: Sets the full intensity matrix with integer bins

    @param data: Raw GCMS data
    @type bin_size: pyms.GCMS.GCMS_data

    @return: Binned IntensityMatrix object
    @rtype: pyms.GCMS.IntensityMatrix

    @author: Qiao Wang
    @author: Andrew Isaac
    @author: Vladimir Likic
    """

    if not isinstance(data, GCMS_data):
        error("data must be an GCMS_data object")

    min_mass = int(round(data.get_min_mass()))
    max_mass = int(round(data.get_max_mass()))
    return __fill_bins(data, min_mass, max_mass, 1)

def __fill_bins(data, min_mass, max_mass, bin_size):

    """
    @summary: Fills the intensity values for all bins

    @param data: Raw GCMS data
    @type bin_size: pyms.GCMS.GCMS_data
    @param min_mass: minimum mass value
    @type min_mass: IntType or FloatType
    @param max_mass: maximum mass value
    @type max_mass: IntType or FloatType
    @param bin_size: bin size
    @type bin_size: IntType or FloatType

    @return: Binned IntensityMatrix object
    @rtype: pyms.GCMS.IntensityMatrix

    @author: Qiao Wang
    @author: Andrew Isaac
    @author: Vladimir Likic
    """

    if not isinstance(data, GCMS_data):
        error("data must be an GCMS_data object")
    if not is_number(min_mass):
        error("'min_mass' must be a number")
    if not is_number(max_mass):
        error("'max_mass' must be a number")
    if not is_number(bin_size):
        error("'bin_size' must be a number")

    # calculate how many bins based on the mass range and bin size
    # round to ensure max mass is included correctly
    num_bins = int(round(float(max_mass - min_mass) / bin_size)) + 1

    #print "Num of bins is:", num_bins

    # initialise masses to bin centres
    mass_list = [i * bin_size + min_mass for i in range(num_bins)]

    # fill the bins
    intensity_matrix = []
    for scan in data.get_scan_list():
        intensitylist = [0.0] * num_bins
        masses = scan.get_mass_list()
        intensities = scan.get_intensity_list()
        for i in range(len(scan)):
            index = int(round(float(masses[i] - min_mass) / bin_size))
            intensitylist[index] += intensities[i]
        intensity_matrix.append(intensitylist)

    return IntensityMatrix(data.get_time_list(), mass_list, intensity_matrix)

def diff(data1, data2):

    """
    @summary: Check if 2 data sets are the same

    @param im1: GCMS data set 1
    @type im1: pyms.GCMS.Class.GCMS_data
    @param im2: GCMS data set 2
    @type im2: pyms.GCMS.Class.GCMS_data

    @author: Qiao Wang
    @author: Andrew Isaac
    @author: Vladimir Likic
    """

    # calculate the time RMSD
    time1 = data1.get_time_list()
    time2 = data2.get_time_list()
    len_time1 = str(len(time1))
    len_time2 = str(len(time2))

    if not len(time1)==len(time2):
        error_type = "Error: the number of time points are different. "
        error_content = "len(time1), and len(time2) are: "
        error(error_type + error_content + len_time1 + ", " + len_time2)
    sum = 0.0
    for i in range(len(time1)):
        sum = sum + (time1[i] - time2[i]) ** 2
    time_RMSD = math.sqrt(sum / len(time1))

    print "Time RMSD: ", time_RMSD
    if time_RMSD > 2.1e-3:
        error("Error: RMSD for time must not exceed the 2.1e - 3")

    #calculate the mass RMSD
    scan_list1 = data1.get_scan_list()
    scan_list2 = data2.get_scan_list()
    len_scan1 = str(len(scan_list1))
    len_scan2 = str(len(scan_list2))

    if not len(scan_list1)==len(scan_list2):
        error_type = "Error: the number of scans are different. "
        error_content = "len(scan_list1), and len(scan_list2) are: "
        error(error_type + error_content + len_scan1 + ", " + len_scan2)

    for j in range(len(scan_list1)):
        mass1 = scan_list1[j].get_mass_list()
        mass2 = scan_list2[j].get_mass_list()
        len_mass1 = str(len(mass1))
        len_mass2 = str(len(mass2))
        if not len(mass1)==len(mass2):
            error_type = "Error: the number of masses are different. "
            error_content = "Scan number: " + str(j) + ". " + "len(mass1), len(mass2) are: "
            error(error_type + error_content + len_mass1 + ", " + len_mass2)
        else:
            mass_RMSD = 0.0
            sum = 0.0
            for k in range(len(mass1)):
                sum = sum + (mass1[k] - mass2[k]) ** 2
            mass_RMSD = math.sqrt(sum / len(mass1))
            print "Mass RMSD for scan point " + str(j) + ": ", mass_RMSD
            if mass_RMSD > 2.0e-5:
                error("Scan number: " + str(j) + ". Error: RMSD for mass must not exceed the 2.0e - 5")
            else:
                pass
        mass1 = []
        mass2 = []

    print "Two data set are the same"