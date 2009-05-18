"""
BillerBiemann deconvolution algorithm
Provides a class to perform Biller and Biemann deconvolution
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

from pyms.Utils.Error import error
from pyms.Utils.Utils import is_list, is_number, is_int
from pyms.GCMS.Class import IonChromatogram

#######################
# structure
# 1) find local maxima per ion, store intensity and scan index
# 2) sum across N scans to compensate for scan type
# 3) sum ions belonging to each maxima scan
#######################

def BillerBiemann(im, window=1):

    """
    @summary: BillerBiemann Deconvolution

        Deconvolution based on the algorithm of Biller and Biemann (1974)

    @param ic: An IonChromatogram object
    @type ic: pyms.Peak.List.Class.Alignment
    @param noise: A number. The noise estimate
    @type noise: FloatType
    @param window: A string. The width of the window for peak
        detection specified as <NUMBER>s or <NUMBER>m, giving a time
        in seconds or minutes, respectively
    @type window: StringType
    @param peak_factor: A number. The value that is multipled to the
        noise level to determine the peak threshold for minimum peak
        intensities
    @type peak_factor: IntType

    @return: A list.
    @rtype: ListType

    @author: Andrew Isaac
    """

    maxima_im = get_maxima_matrix(im)
    sums = []
    numrows = len(maxima_im)
    half = int(window/2)
    for row in range(numrows):
        val = 0
        for ii in range(window):
            if row - half + ii >= 0 and row - half + ii < numrows:
                val += maxima_im[row - half + ii].sum()
        sums.append(val)
    tic = IonChromatogram(numpy.array(sums), im.get_time_list())

    return tic

def get_maxima_indices(ion_intensities):

#    if not is_list(ion_intensities) or not is_number(ion_intensities[0]):
#        error("'ion_intensities' must be a List of numbers")

##TODO: what is tolerance on floating point comparison?

    # find peak inflection points
    # use a 3 point window
    # for a plateau after a rise, need to check if it is the left edge of
    # a peak
    peak_point = []
    edge = -1
    for index in range(len(ion_intensities)-2):
        left = ion_intensities[index]
        mid = ion_intensities[index+1]
        right = ion_intensities[index+2]
        # max in middle
        if mid > left and mid > right:
            peak_point.append(index+1)
            edge = -1  # ignore previous rising edge
        # flat from rise (left of peak?)
        if mid > left and mid == right:
            edge = index+1  # ignore previous rising edge, update latest
        # fall from flat
        if mid == left and mid > right:
            if edge > -1:
                centre = (edge+index+1)/2  # mid point
                peak_point.append(centre)
            edge = -1

    return peak_point

def get_maxima_list(ic):
    peak_point = get_maxima_indices(ic.get_intensity_array())
    mlist = []
    for index in range(len(peak_point)):
        rt = ic.get_time_at_index(peak_point[index])
        intens = ic.get_intensity_at_index(peak_point[index])
        mlist.append([rt, intens])
    return mlist

def get_maxima_matrix(im):
    numrows, numcols = im.get_size()
    # zeroed matrix, size numrows*numcols
    # bad initialisation!! makes numrows copy of the same row!!!
    # maxima_im = [[0]*numcols]*numrows
    maxima_im = numpy.zeros((numrows, numcols))
    raw_im = numpy.array(im.get_matrix_list())

    for col in range(numcols):  # assume all rows have same width
        # 1st, find maxima
        maxima = get_maxima_indices(raw_im[:,col])
        # 2nd, fill intensities
        for row in maxima:
            maxima_im[row, col] = raw_im[row, col]

    return maxima_im
