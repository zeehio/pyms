"""
Provides a class to model signal peak
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

from pyms.GCMS.Class import IntensityMatrix
from pyms.Utils.Error import error
from pyms.Utils.Utils import is_int, is_number, is_list, is_boolean, is_str
from pyms.Utils.IO import open_for_writing, close_for_writing

class Peak:

    """
    @summary: Models a signal peak

    A signal peak object
    A peak object is initialised with retention time and raw peak area.

    @author: Vladimir Likic
    @author: Andrew Isaac
    """

    def __init__(self, rt=0.0, raw_area=0, minutes=False):

        """
        @param rt: Retention time
        @type rt: A number
        @param raw_area: Raw peak area
        @type raw_area: A number
        @param minutes: Retention time units flag. If True, retention time
            is in minutes; if Flase retention time is in seconds
        @type minutes: BooleanType
        """

        if not is_number(rt):
            error("'rt' must be a number")

        if not is_number(raw_area):
            error("'raw_area' must be a number")

        if minutes:
            rt = rt*60.0

        # basic peak attributes
        self.rt = float(rt)
        self.raw_area = float(raw_area)
        self.norm_area = None

        self.tag = None

##TODO: should attributes be private?
        # single ion chromatogram properties
        self.__ic_mass = None

        # these three attributes are required for
        # setting the peak mass spectrum
        self.pt_bounds = None
        self.mass_spectrum = None
        self.mass_list = None

    def set_pt_bounds(self, pt_bounds):

        """
        @summary: Sets peak boundaries in points

        @param pt_bounds: A list containing left, apex, and right
            peak boundaries in points
        @type pt_bounds: ListType

        @return: none
        @rtype: NoneType
        """

        if not is_list(pt_bounds):
            error("'pt_bounds' must be a list")

        if not len(pt_bounds) == 3:
            error("'pt_bounds' must have exactly 3 elements")
        else:
            for item in pt_bounds:
                if not is_int(item):
                    error("'pt_bounds' element not an integer")

        self.pt_bounds = pt_bounds

    def set_ic_mass(self, mz=None):

        """
        @summary: Sets the mass for a single ion chromatogram peak
            Clears the mass spectrum

        @param mz: The mass of the ion chromatogram that the peak is from
        @type mz: FloatType

        @return: none
        @rtype: NoneType
        """

        if not is_number(mz):
            error("'mz' must be a number")
        self.__ic_mass = mz
        # clear mass spectrum
        self.mass_spectrum = None

    def get_ic_mass(self):

        """
        @summary: Gets the mass for a single ion chromatogram peak

        @return: The mass of the single ion chromatogram that the peak is from
        @rtype: FloatType or IntType
        """

        return self.__ic_mass

    def set_mass_spectrum(self, data, from_bounds=False):

        """
        @summary: Sets peak mass spectrum
            Clears the single ion chromatogram mass

        @param data: An IntensityMatrix object
        @type data: pyms.GCMS.IntensityMatrix
        @param from_bounds: Indicator whether to use the attribute
            'pt_bounds' or to find the peak apex from the peak
            retention time
        @type from_bounds: BooleanType

        @return: none
        @rtype: NoneType
        """

        if not isinstance(data, IntensityMatrix):
            error("'data' must be an IntensityMatrix")

        if not is_boolean(from_bounds):
            error("'from_bounds' must be boolean")

        if from_bounds:
            if self.pt_bounds == None:
                error("pt_bounds not set for this peak")
            else:
                pt_apex = self.pt_bounds[1]
        else:
            # get the index of peak apex from peak retention time
            pt_apex = data.get_index_at_time(self.rt)

        # set the mass spectrum
        self.mass_spectrum = data.get_scan_at_index(pt_apex)

        # set mass list for this peak
        self.mass_list = data.get_mass_list()

        # 'mass_spectrum' and 'mass_list' must be consistent.
        if not (len(self.mass_spectrum) == len(self.mass_list)):
            error("mass spectrum data inconsistent")
        # clear single ion chromatogram mass
        self.__ic_mass = None

    def crop_mass_spectrum(self, mass_min, mass_max, verbose=False):

        """
        @summary: Crops mass spectrum

        @param mass_min: Minimum mass value
        @type mass_min: IntType or FloatType
        @param mass_max: Maximum mass value
        @type mass_max: IntType or FloatType
        @param verbose: A verbose flag
        @type verbose: BooleanType

        @return: none
        @rtype: NoneType
        """

        if self.mass_spectrum == None:
            error("mass spectrum not set for this peak")

        if not is_number(mass_min) or not is_number(mass_max):
            error("'mass_min' and 'mass_max' must be numbers")

        new_mass_list = []
        new_mass_spectrum = []

        for ii in range(len(self.mass_list)):

            mass = self.mass_list[ii]
            intensity =  self.mass_spectrum[ii]

            if mass >= mass_min and mass <= mass_max:
                new_mass_list.append(mass)
                new_mass_spectrum.append(intensity)

        self.mass_list = new_mass_list
        self.mass_spectrum = new_mass_spectrum

        if len(new_mass_list) == 0:
            error("mass spectrum empty")
        elif len(new_mass_list) < 10:
            print " WARNING: peak mass spectrum contains < 10 points"

        if verbose:
            print " [ Mass spectrum cropped to %d points:" % \
                    (len(new_mass_list)),
            print "masses %d to %d ]" % (min(new_mass_list),
                    max(new_mass_list))

    def set_peak_tag(self, tag):

        """
        @summary: Sets the peak tag

        @param tag: The supplied peak tag must be either "BLANK" or RF-type
            tag for a reference peak (for example, "RF-SI")
        @type tag: StringType

        @return: none
        @rtype: NoneType
        """

        peak_tag = None

        if tag != None:
            if tag[:3] == "rf-" and tag[3:] != "":
                peak_tag = tag
            elif tag == "blank":
                peak_tag = tag
            else:
                error("incorrect reference peak tag '%s'" % (tag))

        self.tag = peak_tag

    def write_mass_spectrum(self, file_name):

        """
        @summary: Writes the peak mass spectrum to a file

        @param file_name: File for writing the mass spectrum
        @type file_name: StringType

        @return: none
        @rtype: NoneType
        """

        if self.mass_spectrum == None:
            error("mass spectrum for this peak not set")

        if not is_str(file_name):
            error("'file_name' must be a string")

        fp = open_for_writing(file_name)

        for ix in range(len(self.mass_list)):

            fp.write("%8.3f %16.3f\n" % (self.mass_list[ix], \
                    self.mass_spectrum[ix]))

        close_for_writing(fp)

