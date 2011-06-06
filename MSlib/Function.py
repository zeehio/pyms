"""
Functions to create mass spectral library
"""

 #############################################################################
 #                                                                           #
 #    PyMS software for processing of metabolomic mass-spectrometry data     #
 #    Copyright (C) 2005-2011 Vladimir Likic                                 #
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

import sys
sys.path.append('/x/PyMS/')

import string

from pyms.MSlib.Class import msrecord, mslib
from pyms.Utils.IO import file_lines, dump_object, load_object

def load_nist(file_name):    
        
    """
    @summary: Parses through the records of the NIST file
    and builds a mass spectral library

    @param file_name: Name of the NIST file
    @type file_name: StringType

    @return: Mass spectral library
    @rtype: pyms.MSlib.Class.mslib
    
    @author: Saravanan Dayalan
    """
        
    __CMPD_NAME_KEYWORD = "##CAS NAME"
    __CAS_REG_KEYWORD = "##CAS REGISTRY NO"
    __NUM_PEAKS_KEYWORD = "##NPOINTS"
    __RECORD_START_KEYWORD = "##TITLE"

    nist_lines = file_lines(file_name)

    collect_ms_flag = False
    cmpd_counter = 0
    ms_lib = mslib()

    for line in nist_lines:

        fields = string.split(line, "=")

        if len(fields)>0 and fields[0].upper() == __RECORD_START_KEYWORD:
            collect_ms_flag = False

        elif len(fields)>0 and fields[0].upper() == __CMPD_NAME_KEYWORD:
            ms_record = msrecord()
            cmpd_counter = cmpd_counter+1

            if cmpd_counter>1:
                ms_record.set(cmpd_name, cmpd_regno, mi_dict)
                ms_lib.addrecord(ms_record)

            keyword_value = fields[1]
            cmpd_name = keyword_value.strip()

        elif len(fields)>0 and fields[0].upper() == __CAS_REG_KEYWORD:
            keyword_value = fields[1]
            cmpd_regno = keyword_value.strip()

        elif len(fields)>0 and fields[0].upper() == __NUM_PEAKS_KEYWORD:
            keyword_value = fields[1]
            num_peaks_str = keyword_value.strip()
            num_peaks = int(num_peaks_str)
            collect_ms_flag = True
            mi_dict = {}

        elif collect_ms_flag:
            if fields[0][0]!='#':
                mi = string.split(line)
                mass = int(mi[0])
                intensity = int(mi[1])
                mi_dict[mass] = intensity

    ms_record.set(cmpd_name, cmpd_regno, mi_dict)
    ms_lib.addrecord(ms_record)

    return ms_lib


def write_ms_lib(ms_lib, filename):
        
    """
    @summary: Writes the mass spectral library as
    an object to a file

    @param ms_lib: An instance of a mass spectral library
    @type ms_lib: pyms.MSlib.Class.mslib

    @param filename: Name of the file to write the mass
    spectral library as an object
    @type filename: StringType
    """
    
    dump_object(ms_lib, filename)


def load_ms_lib(filename):
        
    """
    @summary: Loads the mass spectral library from the object
    file

    @param filename: Name of the file that stores the mass
    spectral library as an object
    @type filename: StringType

    @return: Mass spectral library
    @rtype: pyms.MSlib.Class.mslib
    """
    
    print "Loading ms library object..."
    lib = load_object(filename)
    return lib


def read_ms_lib(ms_lib):
        
    """
    @summary: Reads the mass spectral library and prints
    the sequence number and the name of each compound 

    @param ms_lib: An instance of a mass spectral library
    @type ms_lib: pyms.MSlib.Class.mslib

    """
        
    counter = 0
    for x in ms_lib.msrecord_list:
        counter = counter + 1
        print counter, x.name    

