"""
Classes to model GC-MS data
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

import numpy
import math
import copy
import operator

from pyms.Utils.Error import error
from pyms.Utils.Utils import is_str, is_int, is_array, is_list, is_number
from pyms.Utils.IO import open_for_writing, close_for_writing
from pyms.Utils.Math import mean, std, median
from pyms.Utils.Time import time_str_secs 

class GCMS_data(object):

    """
    @summary: Generic object for GC-MS data. Contains raw data
        as a list of scans and times

    @author: Qiao Wang
    @author: Andrew Isaac
    @author: Vladimir Likic
    """

    def __init__(self, time_list, scan_list):

        """
        @summary: Initialize the GC-MS data

        @param time_list: List of scan retention times
        @type time_list: ListType
        @param scan_list: List of Scan objects
        @type scan_list: ListType

        @author: Qiao Wang
        @author: Andrew Isaac
        @author: Vladimir Likic
        """

        if not is_list(time_list) or not is_number(time_list[0]):
            error("'time_list' must be a list of numbers")

        if not is_list(scan_list) or not isinstance(scan_list[0], Scan):
            error("'scan_list' must be a list of Scan objects")

        self.__set_time(time_list)
        self.__scan_list = scan_list
        self.__set_min_max_mass()
        self.__calc_tic()

    def __set_time(self, time_list):

        """
        @summary: Sets time-related properties of the data

        @param time_list: List of retention times
        @type time_list: ListType

        @author: Vladimir Likic
        """

        # calculate the time step, its spreak, and along the way
        # check that retention times are increasing
        time_diff_list = []

        for ii in range(len(time_list)-1):
            t1 = time_list[ii] 
            t2 = time_list[ii+1] 
            if not t2 > t1:
                error("problem with retention times detected")
            time_diff = t2 - t1
            time_diff_list.append(time_diff)

        time_step = mean(time_diff_list)
        time_step_std = std(time_diff_list)

        self.__time_list = time_list
        self.__time_step = time_step
        self.__time_step_std = time_step_std
        self.__min_rt = min(time_list)
        self.__max_rt = max(time_list)

    def __set_min_max_mass(self):

        """
        @summary: Sets the min and max mass value

        @author: Qiao Wang
        @author: Andrew Isaac
        @author: Vladimir Likic
        """

        mini = self.__scan_list[0].get_min_mass()
        maxi = self.__scan_list[0].get_max_mass()
        for scan in self.__scan_list:
            tmp_mini = scan.get_min_mass()
            tmp_maxi = scan.get_max_mass()
            if tmp_mini < mini:
                mini = tmp_mini
            if tmp_maxi > maxi:
                maxi = tmp_maxi
        self.__min_mass = mini
        self.__max_mass = maxi

    def get_min_mass(self):

        """
        @summary: Get the min mass value over all scans

        @return: The minimum mass of all the data
        @rtype: FloatType

        @author: Qiao Wang
        @author: Andrew Isaac
        @author: Vladimir Likic
        """

        return self.__min_mass

    def get_max_mass(self):

        """
        @summary: Get the max mass value over all scans

        @return: The maximum mass of all the data
        @rtype: FloatType

        @author: Qiao Wang
        @author: Andrew Isaac
        @author: Vladimir Likic
        """

        return self.__max_mass

    def get_index_at_time(self, time):

        """
        @summary: Returns the nearest index corresponding to the given time

        @param time: Time in seconds
        @type time: FloatType

        @return: Nearest index corresponding to given time
        @rtype: IntType

        @author: Lewis Lee
        @author: Tim Erwin
        @author: Vladimir Likic
        """

        if not is_number(time):
            error("'time' must be a number")

        if time < self.__min_rt or time > self.__max_rt:
            error("time %.2f is out of bounds (min: %.2f, max: %.2f)" %
                  (time, self.__min_rt, self.__max_rt))

        time_list = self.__time_list
        time_diff_min = self.__max_rt
        ix_match = None

        for ix in range(len(time_list)):

            time_diff = math.fabs(time-time_list[ix])

            if time_diff < time_diff_min:
                ix_match = ix
                time_diff_min = time_diff

        return ix_match

    def get_time_list(self):

        """
        @summary: Returns the list of each scan retention time

        @return: A list of each scan retention time
        @rtype: ListType

        @author: Qiao Wang
        @author: Andrew Isaac
        @author: Vladimir Likic
        """

        return copy.deepcopy(self.__time_list)

    def get_scan_list(self):

        """
        @summary: Return a list of the scan objects

        @return: A list of scan objects
        @rtype: ListType

        @author: Qiao Wang
        @author: Andrew Isaac
        @author: Vladimir Likic
        """

        return copy.deepcopy(self.__scan_list)

    def get_tic(self):

        """
        @summary: Returns the total ion chromatogram

        @return: Total ion chromatogram
        @rtype: pyms.GCMS.IonChromatogram

        @author: Andrew Isaac
        """

        return self.__tic

    def __calc_tic(self):
        """
        @summary: Calculate the total ion chromatogram

        @return: Total ion chromatogram
        @rtype: pyms.GCMS.IonChromatogram

        @author: Qiao Wang
        @author: Andrew Isaac
        @author: Vladimir Likic
        """

        intensities = []
        for scan in self.__scan_list:
            intensities.append(sum(scan.get_intensity_list()))
        ia = numpy.array(intensities)
        rt = copy.deepcopy(self.__time_list)
        tic = IonChromatogram(ia, rt)

        self.__tic = tic

    def __scan_values(self, N, first_scan, last_scan):


        """
        @summary: Filters first and last scan values

        @param N: Total number of scans
        @type N: IntType
        @param first_scan: The first scan number
        @type first_scan: IntType
        @param last_scan: The first scan number
        @type last_scan: IntType

        @return: A tuple first_scan, last_scan checked
        @rtype: TupleType

        @author: Vladimir Likic
        """

        if first_scan == None:
           first_scan = 0
        # scan is given. first check that this is sensible
        elif not is_int(first_scan) or first_scan < 1:
           error("first scan must be an integer > 0")
        # adjust for python indexing
        else:
           first_scan = first_scan - 1

        if last_scan == None:
           last_scan = N-1
        # scan is given. first check that this is sensible
        elif not is_int(last_scan) or last_scan < 1:
           error("last scan must be an integer > 0")
        # check that last scan is not greater than the number of scans,
        # and if ok adjust for python indexing
        else:
           if last_scan > N:
               error("last scan greater than max number of scans (%d)" % (N))
           else:
               last_scan = last_scan - 1

        # last_scan can be equal to first_scan, but not greater
        if last_scan < first_scan:
           error("last scan cannot be greater than first scan")

        return first_scan, last_scan

    def trim(self, begin=None, end=None):

        """
        @summary: trims data in the time domain 

        @param begin: begin parameter designating start time or
            scan number 
        @type begin: IntType or StrType
        @param end: end parameter designating start time or
            scan number 
        @type end: IntType or StrType

        The arguments 'begin' and 'end' can be either integers
        (in which case they are taken as the first/last scan
        number for trimming) or strings in which case they are
        treated as time strings and converted to scan numbers.

        At present both 'begin' and 'end' must be of the same
        type, either both scan numbers or time strings.

        @author: Vladimir Likic
        """

        # trim called with defaults, or silly arguments
        if begin == None and end == None:
            print "Nothing to do."
            return # exit immediately
        # trim called with integer arguments
        elif is_int(begin) and is_int(end):
            first_scan = begin
            last_scan = end
        # trim called with time strings as arguments
        elif is_str(begin) and is_str(end):
            start_rt = time_str_secs(begin)
            end_rt = time_str_secs(end)
            # get index and pretend that scan numbers counting
            # from 1 were given, to be consistent with manually
            # given scan numbers
            first_scan = self.get_index_at_time(start_rt) + 1
            last_scan = self.get_index_at_time(end_rt) + 1
        else:
            error("invalid 'begin' and 'end' arguments")

        # get the first and last scans for trimming
        N = len(self.__scan_list)
        # sanity check and reduce scan numbers to python indices
        first_scan, last_scan = self.__scan_values(N, first_scan, last_scan)

        print "Trimming data to between %d and %d scans" % \
                (first_scan+1, last_scan+1)

        scan_list_new = []
        time_list_new = []
        for ii in range(len(self.__scan_list)):
            if ii >= first_scan and ii <= last_scan:
                scan = self.__scan_list[ii]
                time = self.__time_list[ii]
                scan_list_new.append(scan)
                time_list_new.append(time)

        self.__time_list = time_list_new
        self.__scan_list = scan_list_new

        self.__min_rt = min(time_list_new)
        self.__max_rt = max(time_list_new)

    def info(self, print_scan_n=False):

        """
        @summary: Prints some information about the data

        @param print_scan_n: If set to True will print the number
            of m/z values in each scan
        @type print_scan_n: BooleanType

        @author: Vladimir Likic
        """

        # print the summary of simply attributes
        print " Data retention time range: %.3f min -- %.3f min" % \
                (self.__min_rt/60.0, self.__max_rt/60)
        print " Time step: %.3f s (std=%.3f s)" % ( self.__time_step, \
                self.__time_step_std )
        print " Number of scans: %d" % ( len(self.__scan_list) )
        print " Minimum m/z measured: %.3f" % ( self.__min_mass )
        print " Maximum m/z measured: %.3f" % ( self.__max_mass )

        # calculate median number of m/z values measured per scan
        n_list = []
        for ii in range(len(self.__scan_list)):
            scan = self.__scan_list[ii]
            n = len(scan)
            n_list.append(n)
            if print_scan_n: print n
        mz_mean = mean(n_list)
        mz_median = median(n_list)
        print " Mean number of m/z values per scan: %d" % ( mz_mean )
        print " Median number of m/z values per scan: %d" % ( mz_median )

    def write(self, file_root):

        """
        @summary: Writes the entire raw data to two files, one
            'file_root'.I (intensities) and 'file_root'.mz (m/z
            values.

        @param file_root: The rood for the output file names
        @type file_root: StringType
        
        @author: Vladimir Likic
        """

        if not is_str(file_root):
            error("'file_root' must be a string")

        file1 = file_root + ".I"
        file2 = file_root + ".mz"

        fp1 = open_for_writing(file1)
        fp2 = open_for_writing(file2)

        for ii in range(len(self.__scan_list)):

            scan = self.__scan_list[ii]

            intensities = scan.get_intensity_list()
            masses = scan.get_mass_list()

            for item in intensities:
                fp1.write("%8.4f " % (item))
            fp1.write("\n")

            for item in masses:
                fp2.write("%8.4f " % (item))
            fp2.write("\n")

        close_for_writing(fp1)
        close_for_writing(fp2)

    def write_intensities(self, file_name, begin=None, end=None):

        """
        @summary: Writes all intensities to a file

        @param file_name: Output file name
        @type file_name: StringType
        @param begin: The first scan to write
        @type begin: IntType (default: None)
        @param end: The last scan to write
        @type end: IntType (default: None)

        This function loop over all scans, and for each
        scan writes intensities to the file, one intenisity
        per line. Intensities from different scans are
        joined without any delimiters.

        Parameters 'begin' and 'end' are scan numbers not
        python indices, ie. begin = 1 refers to the first
        scan.

        @author: Vladimir Likic
        """

        if not is_str(file_name):
            error("'file_name' must be a string")

        N = len(self.__scan_list)
        first_scan, last_scan = self.__scan_values(N, begin, end)

        # checking arguments passed
        print(" -> Writing scans %d to %d to a file" % (first_scan+1, last_scan+1))

        fp = open_for_writing(file_name)

        for ii in range(len(self.__scan_list)):
            if ii >= first_scan and ii <= last_scan:
                scan = self.__scan_list[ii]
                intensities = scan.get_intensity_list()
                for I in intensities:
                    fp.write("%8.4f\n" % (I))

        close_for_writing(fp)

class Scan(object):

    """
    @summary: Generic object for a single Scan's raw data

    @author: Qiao Wang
    @author: Andrew Isaac
    @author: Vladimir Likic
    """

    def __init__(self, mass_list, intensity_list):

        """
        @summary: Initialize the Scan data

        @param mass_list: mass values
        @type mass_list: ListType

        @param intensity_list: intensity values
        @type intensity_list: ListType

        @author: Qiao Wang
        @author: Andrew Isaac
        @author: Vladimir Likic
        """

        if not is_list(mass_list) or not is_number(mass_list[0]):
            error("'mass_list' must be a list of numbers")
        if not is_list(intensity_list) or \
           not is_number(intensity_list[0]):
            error("'intensity_list' must be a list of numbers")

        self.__mass_list = mass_list
        self.__intensity_list = intensity_list
        self.__min_mass = min(mass_list)
        self.__max_mass = max(mass_list)

    def __len__(self):

        """
        @summary: Returns the length of the Scan object

        @return: Length of Scan
        @rtype: IntType

        @author: Andrew Isaac
        """

        return len(self.__mass_list)

    def get_mass_list(self):

        """
        @summary: Returns the masses for the current scan

        @return: the masses
        @rtype: ListType

        @author: Qiao Wang
        @author: Andrew Isaac
        @author: Vladimir Likic
        """

        return copy.deepcopy(self.__mass_list)

    def get_intensity_list(self):

        """
        @summary: Returns the intensities for the current scan

        @return: the intensities
        @rtype: ListType

        @author: Qiao Wang
        @author: Andrew Isaac
        @author: Vladimir Likic
        """

        return copy.deepcopy(self.__intensity_list)

    def get_min_mass(self):

        """
        @summary: Returns the minimum m/z value in the scan

        @return: Minimum m/z
        @rtype: Float

        @author: Andrew Isaac
        """

        return self.__min_mass

    def get_max_mass(self):

        """
        @summary: Returns the maximum m/z value in the scan

        @return: Maximum m/z
        @rtype: Float

        @author: Andrew Isaac
        """

        return self.__max_mass

class IntensityMatrix(object):

    """
    @summary: Intensity matrix of binned raw data

    @author: Andrew Isaac
    """

    def __init__(self, time_list, mass_list, intensity_matrix):

        """
        @summary: Initialize the IntensityMatrix data

        @param time_list: Retention time values
        @type time_list: ListType

        @param mass_list: Binned mass values
        @type mass_list: ListType

        @param intensity_list: Binned intensity values per scan
        @type intensity_list: ListType

        @author: Andrew Isaac
        """

        # sanity check
        if not is_list(time_list) or not is_number(time_list[0]):
            error("'time_list' must be a list of numbers")
        if not is_list(mass_list) or not is_number(mass_list[0]):
            error("'mass_list' must be a list of numbers")
        if not is_list(intensity_matrix) or \
           not is_list(intensity_matrix[0]) or \
           not is_number(intensity_matrix[0][0]):
            error("'intensity_matrix' must be a list, of a list, of numbers")
        if not len(time_list) == len(intensity_matrix):
            error("'time_list' is not the same length as 'intensity_matrix'")
        if not len(mass_list) == len(intensity_matrix[0]):
            error("'mass_list' is not the same size as 'intensity_matrix'"
                " width")

        self.__time_list = time_list
        self.__mass_list = mass_list
        self.__intensity_matrix = intensity_matrix

        self.__min_mass = min(mass_list)
        self.__max_mass = max(mass_list)

    def get_size(self):

        """
        @summary: Gets the size of intensity matrix

        @return: Number of rows and cols
        @rtype: IntType

        @author: Qiao Wang
        @author: Andrew Isaac
        @author: Vladimir Likic
        """

        row = len(self.__intensity_matrix)
        col = len(self.__intensity_matrix[0])

        return row, col

    def get_ic_at_index(self, index):

        """
        @summary: Returns the ion chromatogram at the specified index

        @param index: Index of an ion chromatogram in the intensity data matrix
        @type index: IntType

        @return: Ion chromatogram at given index
        @rtype: IonChromatogram

        @author: Qiao Wang
        @author: Andrew Isaac
        @author: Vladimir Likic
        """

        index = int(index)
        try:
            ia = []
            for i in range(len(self.__intensity_matrix)):
                ia.append(self.__intensity_matrix[i][index])
        except IndexError:
            error("index out of bounds.")

        ia = numpy.array(ia)
        mass = self.get_mass_at_index(index)
        rt = copy.deepcopy(self.__time_list)

        return IonChromatogram(ia, rt, mass)

    def get_ic_at_mass(self, mass = None):

        """
        @summary: Returns the ion chromatogram for the specified mass
            The nearest binned mass to mass is used

        If no mass value is given, the function returns the total
        ion chromatogram.

        @param mass: Mass value of an ion chromatogram
        @type mass: IntType

        @return: Ion chromatogram for given mass
        @rtype: IonChromatogram
        """

        if mass == None:
            return self.get_tic()

        if mass < self.__min_mass or mass > self.__max_mass:
            error("mass is out of range")

        index = self.get_index_of_mass(mass)
        return self.get_ic_at_index(index)

    def get_mass_list(self):

        """
        @summary: Returns a list of the bin masses

        @return: Binned mass list
        @rtype: ListType

        @author: Qiao Wang
        @author: Andrew Isaac
        @author: Vladimir Likic
        """

        return copy.deepcopy(self.__mass_list)

    def get_ms_at_index(self, index):

        """
        @summary: Returns a mass spectrum for a given scan index

        @param index: The index of the scan
        @type index: IntType

        @return: Mass spectrum
        @rtype: pyms.GCMS.MassSpectrum

        @author: Andrew Isaac
        """

        scan = self.get_scan_at_index(index)

        return MassSpectrum(self.__mass_list, scan)

    def get_scan_at_index(self, index):

        """
        @summary: Returns the spectral intensities for scan index

        @param index: The index of the scan
        @type index: IntType

        @return: Intensity values of scan spectra
        @rtype: ListType

        @author: Andrew Isaac
        """

        if index < 0 or index >= len(self.__intensity_matrix):
            error("index out of range")

        return copy.deepcopy(self.__intensity_matrix[index])

    def get_min_mass(self):

        """
        @summary: Returns the maximum binned mass

        @return: The maximum binned mass
        @rtype: FloatType

        @author: Andrew Isaac
        """

        return self.__min_mass

    def get_max_mass(self):

        """
        @summary: Returns the maximum binned mass

        @return: The maximum binned mass
        @rtype: FloatType

        @author: Andrew Isaac
        """

        return self.__max_mass

    def get_mass_at_index(self, index):

        """
        @summary: Returns binned mass at index

        @param index: Index of binned mass
        @type index: IntType

        @return: Binned mass
        @rtype: IntType

        @author: Andrew Isaac
        """

        if index < 0 or index >= len(self.__mass_list):
            error("index out of range")

        return self.__mass_list[index]

    def get_index_of_mass(self, mass):

        """
        @summary: Returns the index of mass in the list of masses
            The nearest binned mass to mass is used

        @param mass: Mass to lookup in list of masses
        @type mass: FloatType

        @return: Index of mass closest to given mass
        @rtype: IntType

        @author: Andrew Isaac
        """

        best = self.__max_mass
        index = 0
        for i in range(len(self.__mass_list)):
            tmp = abs(self.__mass_list[i] - mass)
            if tmp < best:
                best = tmp
                index = i
        return index

    def get_matrix_list(self):

        """
        @summary: Returns a copy of the intensity matrix as a
            list of lists of floats

        @return: Matrix of intensity values
        @rtype: ListType

        @author: Andrew Isaac
        """

        return copy.deepcopy(self.__intensity_matrix)

    def get_time_list(self):

        """
        @summary: Returns a copy of the time list

        @return: List of retention times
        @rtype: ListType

        @author: Andrew Isaac
        """

        return copy.deepcopy(self.__time_list)

    def get_index_at_time(self, time):

        """
        @summary: Returns the nearest index corresponding to the given time

        @param time: Time in seconds
        @type time: FloatType

        @return: Nearest index corresponding to given time
        @rtype: IntType

        @author: Lewis Lee
        @author: Tim Erwin
        @author: Vladimir Likic
        """

        if not is_number(time):
            error("'time' must be a number")

        if time < min(self.__time_list) or time > max(self.__time_list):
            error("time %.2f is out of bounds (min: %.2f, max: %.2f)" %
                  (time, self.__min_rt, self.__max_rt))

        time_list = self.__time_list
        time_diff_min = max(self.__time_list)
        ix_match = None

        for ix in range(len(time_list)):

            time_diff = math.fabs(time-time_list[ix])

            if time_diff < time_diff_min:
                ix_match = ix
                time_diff_min = time_diff

        return ix_match

    def reduce_mass_spectra(self, N=5, min_mz_diff=1.5):

        """
        @summary: Reduces mass spectra by retaining top N
        intensities, and discarding all other intensities.

        @param N: The number of top intensities to keep
        @type N: IntType
        @param min_mz_diff: The minimum difference in m/z
            required to keep the intensity. For example,
            if two m/z values with greater intensity are
            77.0 and 78.0, and min_mz_diff=1.5, the intensity
            for for m/z=77 will be discarded.
        @type min_mz_diff: FloatType

        @return: Does not return a value, changes the variable
            self.__intensity_matrix

        @author: Vladimir Likic
        """

        # loop over all mass spectral scans
        for ii in range(len(self.__intensity_matrix)):

            # get next mass spectrum
            intensity_list = self.__intensity_matrix[ii]

            # sort mass spectrum intensities
            tmp_list = sorted(enumerate(intensity_list), key=operator.itemgetter(1))
            tmp_list.reverse()

            # a sorted list of indices and associated dictionary of intensities
            index_list = []
            intensity_dict = {}
            mass_dict = {}
            for item in tmp_list:
                ix = item[0]
                index_list.append(ix)
                intensity_dict[ix] = intensity_list[ix]
                mass_dict[ix] = self.__mass_list[ix]

            # limit to first N intensities
            top_index_list = []            

            cntr = 0

            while len(top_index_list) < N:

                ix_crnt = index_list[cntr]
                ix_crnt_mass  = mass_dict[ix_crnt]

                mass_ok = True

                for ix in top_index_list:
                   ix_mass = mass_dict[ix] 
                   if math.fabs(ix_crnt_mass - ix_mass) < min_mz_diff:
                       mass_ok = False
 
                if mass_ok:
                    top_index_list.append(ix_crnt)

                cntr = cntr + 1

            # build a cleaned up mass spectrum
            N = len(intensity_list)
            intensity_list_new = []

            for jj in range(N):
                if jj in top_index_list:
                    intensity_list_new.append(intensity_dict[jj])
                else:
                    intensity_list_new.append(0.0)

            self.__intensity_matrix[ii] = intensity_list_new
 
## get_ms_at_time()

class IonChromatogram(object):

    """
    @summary: Models ion chromatogram

    An ion chromatogram is a set of intensities as a function of retention
    time. This can can be either m/z channel intensities (for example, ion
    chromatograms at m/z=65), or cumulative intensities over all measured
    m/z. In the latter case the ion chromatogram is total ion chromatogram
    (TIC).

    The nature of an IonChromatogram object can be revealed by inspecting
    the value of the attribute '__mass'. This is se to the m/z value of the
    ion chromatogram, or to None for TIC.

    @author: Lewis Lee
    @author: Vladimir Likic
    """

    def __init__(self, ia, time_list, mass=None):

        """
        @param ia: Ion chromatogram intensity values
        @type ia: numpy.array
        @param time_list: A list of ion chromatogram retention times
        @type time_list: ListType
        @param mass: Mass of ion chromatogram (Null if TIC)
        @type mass: IntType

        @author: Lewis Lee
        @author: Vladimir Likic
        @author: Vladimir Likic
        """

        if not isinstance(ia, numpy.ndarray):
            error("'ia' must be a numpy array")

        if not is_list(time_list) or not is_number(time_list[0]):
            error("'time_list' must be a list of numbers")

        if len(ia) != len(time_list):
            error("Intensity array and time list differ in length")

        self.__ia = ia
        self.__time_list = time_list
        self.__mass = mass
        self.__time_step = self.__calc_time_step(time_list)
        self.__min_rt = min(time_list)
        self.__max_rt = max(time_list)

    def __len__(self):

        """
        @summary: Returns the length of the IonChromatogram object

        @return: Length of ion chromatogram
        @rtype: IntType

        @author: Lewis Lee
        @author: Vladimir Likic
        """

        return self.__ia.size

    def get_intensity_at_index(self, ix):

        """
        @summary: Returns intensity at given index

        @param ix: An index
        @type ix: IntType

        @return: Intensity value
        @rtype: FloatType

        @author: Lewis Lee
        @author: Vladimir Likic
        """

        if not is_int(ix):
            error("index not an integer")

        if ix < 0 or ix > self.__ia.size - 1:
            error("index out of bounds")

        return self.__ia[ix]

    def get_intensity_array(self):

        """
        @summary: Returns the entire intensity array

        @return: Intensity array
        @rtype: numpy.ndarray

        @author: Lewis Lee
        @author: Vladimir Likic
        """

        return self.__ia

    def get_time_at_index(self, ix):

        """
        @summary: Returns time at given index

        @param ix: An index
        @type ix: IntType

        @return: Time value
        @rtype: FloatType

        @author: Lewis Lee
        @author: Vladimir Likic
        """

        if not is_int(ix):
            error("index not an integer")

        if ix < 0 or ix > len(self.__time_list) - 1:
            error("index out of bounds")

        return self.__time_list[ix]

    def get_time_list(self):

        """
        @summary: Returns the time list

        @return: Time list
        @rtype: ListType

        @author: Lewis Lee
        @author: Vladimir Likic
        """

        return self.__time_list

    def get_time_step(self):

        """
        @summary: Returns the time step

        @return: Time step
        @rtype: FloatType

        @author: Lewis Lee
        @author: Vladimir Likic
        """

        return self.__time_step

    def __calc_time_step(self, time_list):

        """
        @summary: Calculates the time step

        @param time_list: A list of retention times
        @type time_list: ListType

        @return: Time step value
        @rtype: FloatType

        @author: Lewis Lee
        @author: Vladimir Likic
        """

        td_list = []
        for ii in range(len(time_list)-1):
            td = time_list[ii+1]-time_list[ii]
            td_list.append(td)

        td_array = numpy.array(td_list)
        time_step = td_array.mean()

        return time_step

    def get_index_at_time(self, time):

        """
        @summary: Returns the nearest index corresponding to the given time

        @param time: Time in seconds
        @type time: FloatType

        @return: Nearest index corresponding to given time
        @rtype: IntType

        @author: Lewis Lee
        @author: Tim Erwin
        @author: Vladimir Likic
        @author: Milica Ng
        """

        if not is_number(time):
            error("'time' must be a number")

        if time < self.__min_rt or time > self.__max_rt:
            error("time %.2f is out of bounds (min: %.2f, max: %.2f)" %
                  (time, self.__min_rt, self.__max_rt))

        time_list = self.__time_list
        time_diff_min = self.__max_rt
        ix_match = None

        for ix in range(len(time_list)):


            time_diff = math.fabs(time-time_list[ix])

            if time_diff < time_diff_min:
                ix_match = ix
                time_diff_min = time_diff

        return ix_match

    def is_tic(self):

        """
        @summary: Returns True if the ion chromatogram is a total ion
            chromatogram (TIC), or False otherwise

        @return: A boolean value indicating if the ion chromatogram
            is a total ion chromatogram (True) or not (False)
        @rtype: BooleanType

        @author: Lewis Lee
        @author: Vladimir Likic
        """

        if self.__mass == None:
            return True
        else:
            return False

    def write(self, file_name, minutes=False):

        """
        @summary: Writes the ion chromatogram to the specified file

        @param file_name: Output file name
        @type file_name: StringType
        @param minutes: A boolean value indicating whether to write
            time in minutes
        @type minutes: BooleanType

        @return: none
        @rtype: NoneType

        @author: Lewis Lee
        @author: Vladimir Likic
        """

        if not is_str(file_name):
            error("'file_name' must be a string")

        fp = open_for_writing(file_name)

        time_list = copy.deepcopy(self.__time_list)

        if minutes:
            for ii in range(len(time_list)):
                time_list[ii] = time_list[ii]/60.0

        for ii in range(len(time_list)):
            fp.write("%8.4f %#.6e\n" % (time_list[ii], self.__ia[ii]))

        close_for_writing(fp)

class MassSpectrum(object):

    """
    @summary: Models a binned mass spectrum

    @author: Andrew Isaac
    @author: Qiao Wang
    @author: Vladimir Likic
    """

    def __init__(self, mass_list, intensity_list):

        """
        @summary: Initialise the MassSpectrum

        @para mass_list: List of binned masses
        @type mass_list: ListType
        @para intensity_list: List of binned intensities
        @type intensity_list: ListType

        @author: Andrew Isaac
        @author: Qiao Wang
        @author: Vladimir Likic
        """

        if not is_list(mass_list) or not is_number(mass_list[0]):
            error("'mass_list' must be a list of numbers")
        if not is_list(intensity_list) or \
           not is_number(intensity_list[0]):
            error("'intensity_list' must be a list of numbers")
        if not len(mass_list) == len(intensity_matrix):
            error("'mass_list' is not the same size as 'intensity_matrix'")

        #TODO: should these be public, or accessed through methods???
        self.mass_list = mass_list
        self.mass_spec = intensity_list

    def __len__(self):

        """
        @summary: Length of the MassSpectrum

        @return: Length of the MassSpectrum (Number of bins)
        @rtype: IntType

        @author: Andrew Isaac
        @author: Qiao Wang
        @author: Vladimir Likic
        """

        return len(self.mass_list)

