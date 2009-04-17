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

from pyms.Utils.Error import error
from pyms.Utils.Utils import is_number, is_str, is_array, is_list, is_int
from pyms.Utils.IO import open_for_writing, close_for_writing
from pyms.Utils.IO import save_data

def export_csv(data, root_name):

    """
    @summary: Exports data to the CSV format

    Calling export_csv("NAME") will create NAME.im.csv, NAME.rt.csv,
    and NAME.mz.csv where these are the intensity matrix, retention
    time vector, and m/z vector.

    @param root_name: Root name for the output files
    @type root_name: StringType

    @return: none
    @rtype: NoneType

    @author: Milica Ng
    """

    if not is_str(root_name):
        error("'root_name' is not a string")

    # export 2D matrix of intensities into CSV format
    im = data.get_intensity_matrix()
    save_data(root_name+'.im.csv', im, sep=",")

    # export 1D vector of m/z's, corresponding to rows of
    # the intensity matrix, into CSV format
    mass_list = data.get_mass_list()
    save_data(root_name+'.mz.csv', mass_list, sep=",")

    # export 1D vector of retention times, corresponding to
    # columns of the intensity matrix, into CSV format
    time_list = data.get_time_list()
    save_data(root_name+'.rt.csv', time_list, sep=",")

def export_leco_csv(data, file_name):

    """
    @summary: Exports data in LECO CSV format

    @param data: Intensity matrix
    @type data: numpy.ndarray
    @param file_name: File name
    @type file_name: StringType

    @return: none
    @rtype: NoneType

    @author: Andrew Isaac
    @author: Vladimir Likic
    """

    if not is_str(file_name):
        error("'file_name' is not a string")

    im = data.get_intensity_matrix()
    mass_list = data.get_mass_list()
    time_list = data.get_time_list()

    fp = open_for_writing(file_name)

    # Format is text header with:
    # "Scan","Time",...
    # and the rest is "TIC" or m/z as text, i.e. "50","51"...
    # The following lines are:
    # scan_number,time,value,value,...
    # scan_number is an int, rest seem to be fixed format floats.  The format
    # is 0.000000e+000

    # write header
    fp.write("\"Scan\",\"Time\"")
    for ii in mass_list:
        if is_number(ii):
            fp.write(",\"%d\"" % int(ii))
        else:
            error("mass list datum not a number")
    fp.write("\r\n")

    # write lines
    for ii in range(len(time_list)):
        fp.write("%s,%#.6e" % (ii, time_list[ii]))
        for jj in range(len(im[ii])):
            if is_number(im[ii][jj]):
                fp.write(",%#.6e" % (im[ii][jj]))
            else:
                error("datum not a number")
        fp.write("\r\n")

    close_for_writing(fp)
