"""
Functions for reading manufacturer specific ANDI-MS data files
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

import math, copy

import numpy

from pyms.GCMS.Class import GCMS_data
from pyms.GCMS.Class import Scan
from pyms.Utils.Error import error, stop
from pyms.Utils.Utils import is_str, is_int, is_float, is_number, is_list
from pyms.Utils.Time import time_str_secs
from pycdf import CDF, CDFError

def ANDI_reader(file_name):

    """
    @summary: A reader for ANDI-MS NetCDF files, returns
        a GC-MS data object

    @param file_name: The name of the ANDI-MS file
    @type file_name: StringType

    @author: Qiao Wang
    @author: Andrew Isaac
    @author: Vladimir Likic
    """

    # the keys used to retrieve certain data from the NetCDF file
    __MASS_STRING = "mass_values"
    __INTENSITY_STRING = "intensity_values"
    __TIME_STRING = "scan_acquisition_time"

    if not is_str(file_name):
        error("'file_name' must be a string")
    try:
        file = CDF(file_name)
    except CDFError:
        error("Cannot open file '%s'" % file_name)

    print " -> Reading netCDF file '%s'" % (file_name)

    scan_list = []
    mass = file.var(__MASS_STRING)
    intensity = file.var(__INTENSITY_STRING)
    mass_values = mass.get().tolist()
    mass_list = []
    mass_previous = mass_values[0]
    mass_list.append(mass_previous)
    intensity_values = intensity.get().tolist()
    intensity_list = []
    intensity_previous = intensity_values[0]
    intensity_list.append(intensity_previous)
    if not len(mass_values) == len(intensity_values):
        error("length of mass_list is not equal to length of intensity_list !")
    for i in range(len(mass_values) - 1):
        # assume masses in ascending order until new scan
        if mass_previous <= mass_values[i + 1]:
            mass_list.append(mass_values[i + 1])
            mass_previous = mass_values[i + 1]
            intensity_list.append(intensity_values[i + 1])
            intensity_previous = intensity_values[i + 1]
        # new scan
        else:
            scan_list.append(Scan(mass_list, intensity_list))
            mass_previous = mass_values[i + 1]
            intensity_previous = intensity_values[i + 1]
            mass_list = []
            intensity_list = []
    # store final scan
    scan_list.append(Scan(mass_list, intensity_list))
    time = file.var(__TIME_STRING)
    time_list = time.get().tolist()

    # sanity check
    if not len(time_list) == len(scan_list):
        error("number of time points does not equal the number of scans")

    data = GCMS_data(time_list, scan_list)

    return data

