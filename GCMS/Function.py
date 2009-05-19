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

import math, sys

from pyms.Utils.Error import error
from pyms.Utils.Utils import is_number, is_str, is_array, is_list, is_int
from pyms.GCMS.Class import GCMS_data, IntensityMatrix, IonChromatogram
from pyms.Utils.Time import time_str_secs
from pyms.Utils.Math import rmsd

def build_intensity_matrix(data, bin_size=1):

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
    @summary: Compares two GCMS_data objects

    @param im1: GCMS data set 1
    @type im1: pyms.GCMS.Class.GCMS_data
    @param im2: GCMS data set 2
    @type im2: pyms.GCMS.Class.GCMS_data

    @author: Qiao Wang
    @author: Andrew Isaac
    @author: Vladimir Likic
    """

    # get time attributes
    time_list1 = data1.get_time_list()
    time_list2 = data2.get_time_list()

    #
    # First, check if two data sets have the same number of retention
    # times.
    #
    if not len(time_list1) == len(time_list2):
        print " -> The number of retention time points different."
        print " First data set: %d time points" % (len(time_list1))
        print " Second data set: %d time points" % (len(time_list2))
        print " Data sets are different."
        return
    else:
        time_rmsd = rmsd(time_list1, time_list2)
        print " Data sets have the same number of time points."
        print "   Time RMSD: %.2e" % ( time_rmsd )

    #
    # Second, check if each scan has the same number of m/z intensities
    #

    print " Checking for consistency in scan lengths ...",
    sys.stdout.flush()

    scan_list1 = data1.get_scan_list()
    scan_list2 = data2.get_scan_list()
    if not len(scan_list1) == len(scan_list2):
        # since the number of rention times are the same, this indicated
        # some unexpected problem with data 
        error("inconsistency in data detected")

    N = len(scan_list1)

    for ii in range(N):
        scan1 = scan_list1[ii]
        scan2 = scan_list2[ii]
        mass_list1 = scan1.get_mass_list()
        mass_list2 = scan2.get_mass_list()
        if len(mass_list1) != len(mass_list2):
            print "\n Different number of points detected in scan no. %d" % ( ii )
            print " Data sets are different."
            return

    print "OK"

    #
    # Third, if here, calculate the max RMSD for m/z and intensities
    #

    print " Calculating maximum RMSD for m/z values and intensities ...",
    sys.stdout.flush()

    max_mass_rmsd = 0.0
    max_intensity_rmsd = 0.0

    for ii in range(N):
        scan1 = scan_list1[ii]
        scan2 = scan_list2[ii]
        mass_list1 = scan1.get_mass_list()
        mass_list2 = scan2.get_mass_list()
        intensity_list1 = scan1.get_intensity_list()
        intensity_list2 = scan2.get_intensity_list()
        mass_rmsd = rmsd(mass_list1, mass_list2)
        if mass_rmsd > max_mass_rmsd:
            max_mass_rmsd = mass_rmsd
        intensity_rmsd = rmsd(intensity_list1, intensity_list2)
        if intensity_rmsd > max_intensity_rmsd:
            max_intensity_rmsd = intensity_rmsd

    print "\n   Max m/z RMSD: %.2e" % ( max_mass_rmsd )
    print "   Max intensity RMSD: %.2e" % ( max_intensity_rmsd )

def is_ionchromatogram(arg):

    """
    @summary: Returns True if the argument is a pyms.IO.Class.IonCromatogram
        object, False otherwise

    @param arg: The argument to be evaluated as IonCromatogram object
    @type arg: arbitrary

    @return: A boolean indicator True or False
    @rtype: BooleanType

    @author: Vladimir Likic
    """

    if isinstance(arg,IonChromatogram):
        return True
    else:
        return False

def ic_window_points(ic, window_sele, half_window=False):

    """
    @summary: Converts window selection parameter into points based on
        the time step in an ion chromatogram 

    @param ic: ion chromatogram object relevant for the conversion
    @type ic: pyms.IO.Class.IonChromatogram
    @param window_sele: The window selection parameter. This can be an
        integer or time string. If integer, taken as the number of points.
        If a string, must of the form "<NUMBER>s" or "<NUMBER>m",
        specifying a time in seconds or minutes, respectively
    @type window_sele: IntType or StringType 
    @param half_window: Specifies whether to return half-window
    @type half_window: Booleantype

    @author: Vladimir Likic
    """

    if not is_int(window_sele) and not is_str(window_sele):
        error("'window' must be either an integer or a string")

    if is_int(window_sele):

        if half_window:
            if window_sele % 2 == 0: 
                error("window must be an odd number of points")
            else: 
                points = int(math.floor(window_sele*0.5))
        else:
            points = window_sele
    else:
        time = time_str_secs(window_sele)
        time_step = ic.get_time_step()

        if half_window:
            time = time*0.5

        points = int(math.floor(time/time_step))

    if half_window:
        if points < 1: error("window too small (half window=%d)" % (points))
    else:
        if points < 2: error("window too small (window=%d)" % (points))

    return points
