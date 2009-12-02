"""
Functions for Flux.MassExtract
"""


 #############################################################################
 #                                                                           #
 #    PyMS software for processing of metabolomic mass-spectrometry data     #
 #    Copyright (C) 2005-9 Vladimir Likic                                    #
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


from pyms.Flux.Class import MIDS
from pyms.Noise.SavitzkyGolay import savitzky_golay
from pyms.Baseline.TopHat import tophat

def peak_apex(ic, rt, win_size):

    """
    @summary: Locate peak apex, as a local maximum 

    @param ic: Ion Chromatogram
    @type ic: pyms.GCMS.Class.IonChromatogram
    @param rt: Retention time
    @type rt: FloatType
    @param win_size: Window for a local maximum search
    @type win_size: FloatType
    
    @return: Ion Chromatogram index 
    @rtype: IntType
    
    @author: Milica Ng
    """

    # get indices at retention time +- 1/2 window size
    lb_rt = rt - (win_size/2)
    lb_ix = ic.get_index_at_time(lb_rt)
    ub_rt = rt + (win_size/2)
    ub_ix = ic.get_index_at_time(ub_rt)

    # find peak_apex (local maximum between lb & ub indices)
    peak_apex_int = ic.get_intensity_at_index(lb_ix)
    peak_apex_ix = lb_ix

    for ix in range(lb_ix,ub_ix+1):
       if ic.get_intensity_at_index(ix) > peak_apex_int:
          peak_apex_int = ic.get_intensity_at_index(ix)
          peak_apex_ix = ix
    # peak_apex_rt = ic.get_time_at_index(peak_apex_ix)

    return peak_apex_ix

def left_boundary(ic, peak_apex_ix, peak_apex_int):

    """
    @summary: Locate left boundary, as a local minimum 

    @param ic: Ion Chromatogram
    @type ic: pyms.GCMS.Class.IonChromatogram
    @param peak_apex_ix: Peak apex index
    @type peak_apex_ix: IntType
    @param peak_apex_int: Peak apex intensity
    @type peak_apex_ix: FloatType
    
    @return: Left boundary index 
    @rtype: IntType
    
    @author: Milica Ng
    """

    # find left_boundary (first local minimum to the left of peak apex)
    left_boundary_int = peak_apex_int
    ix = peak_apex_ix
    while ix >= 0:
        if ic.get_intensity_at_index(ix) > left_boundary_int:
            break
        else:
            left_boundary_int = ic.get_intensity_at_index(ix)
            ix = ix - 1
    left_boundary_ix = ix + 1
    # left_boundary_rt = ic.get_time_at_index(left_boundary_ix)   

    return left_boundary_ix

def right_boundary(ic, peak_apex_ix, peak_apex_int):

    """
    @summary: Locate right boundary, as a local minimum 

    @param ic: Ion Chromatogram
    @type ic: pyms.GCMS.Class.IonChromatogram
    @param peak_apex_ix: Peak apex index
    @type peak_apex_ix: IntType
    @param peak_apex_int: Peak apex intensity
    @type peak_apex_ix: FloatType
    
    @return: Right boundary index 
    @rtype: IntType
    
    @author: Milica Ng
    """

    # find right_boundary (first local minimum to the right of peak apex)
    right_boundary_int = peak_apex_int
    ix = peak_apex_ix
    while ix < len(ic):
        if ic.get_intensity_at_index(ix) > right_boundary_int:
            break
        else:
            right_boundary_int = ic.get_intensity_at_index(ix)
            ix = ix + 1
    right_boundary_ix = ix - 1
    # right_boundary_rt = ic.get_time_at_index(right_boundary_ix)
    return right_boundary_ix


def calculate_mid(ic_list, ave_left_ix, ave_right_ix):

    """
    @summary: Calculate mass isotopomer distribution

    @param ion: Ion m/z
    @type ion: IntType
    @param mid_size: mid_size (n+1) is the number of masses (M, M+1, ..., M+n)
    @type mid_size: IntType
    @param ave_left_ix: Average left boundary index for the ions
    @type ave_left_ix: IntType
    @param ave_right_ix: Average right boundary index for the ions
    @type ave_right_ix: IntType
    @param out_file: Name of the output file
    @type out_file: StringType
    
    @return: None
    @rtype: NoneType
    
    @author: Milica Ng
    """

    # calculate mass isotopomer distribution
    mid = []
    for ic in ic_list:
        # get ion chromatogram at specified m/z
        # ic = im.get_ic_at_mass(mz)
        area = 0
        # integrate (sum intensity from average left_boundary_ix to average right_boundary_ix 
        for ix in range(ave_left_ix, ave_right_ix+1):
            area = area + ic.get_intensity_at_index(ix)
        mid.append(area)

    return mid

def extract_mid(file_name, im, mids, win_size, noise):

    """
    @summary: Method for extracting mass isotopomer distribution (MID)

    @param file_name: File number
    @type file_name: StringType
    @param im: Intensity matrix
    @type im: pyms.GCMS.Class.IntensityMatrix
    @param mids: Mass isotopomer distribution data
    @type mids: pyms.Flux.Class.MIDS
    @param mid_size: Size of the mass distribution vector
    @type mid_size: IntType
    @param win_size: Size of the window for local maximum search
    @type win_size: FloatType

    @return: Mass isotopomer distribution data
    @rtype: pyms.Flux.Class.MIDS
    
    @author: Milica Ng
    """
    
    # get name, rt, ion and MID size
    compound = mids.get_name()
    # print compound,
    rt = mids.get_rt()
    # print ' rt:', rt,
    ion = mids.get_ion()
    # print ' ion:', ion,
    mid_size = mids.get_mid_size()
    # print ' mid_size:', mid_size

    # initialise variables
    left_boundary_sum = 0
    right_boundary_sum = 0
    sum_count = 0
    ic_list = []

    # loop over required number of mass isotopomers
    for mz in range(ion, ion+mid_size):

        # get ion chromatogram at current m/z
        # print 'about to get ic at mass', mz
        # todo: check if mass is out of range, raise an error!
        ic = im.get_ic_at_mass(mz) 
        # smooth noise + correct baseline
        ic_smooth = savitzky_golay(ic)
        ic_bc = tophat(ic_smooth, struct="1.5m")
        ic = ic_bc
        # store ic's to pass them to calculate_mid
        ic_list.append(ic)

        # get indices for apex and left/right boundary
        peak_apex_ix = peak_apex(ic, rt, win_size)
        peak_apex_int = ic.get_intensity_at_index(peak_apex_ix) 
        left_boundary_ix = left_boundary(ic, peak_apex_ix, peak_apex_int)
        right_boundary_ix = right_boundary(ic, peak_apex_ix, peak_apex_int)

        # sum left and right boundary indices in the current file
        if peak_apex_int > noise: 
            left_boundary_sum = left_boundary_sum + left_boundary_ix
            right_boundary_sum = right_boundary_sum + right_boundary_ix
            sum_count = sum_count + 1

    if sum_count > 0:

        # find average of left and right boundary index in the current file
        ave_left_ix = left_boundary_sum/sum_count
        ave_right_ix = right_boundary_sum/sum_count

        # calculate mass isotopomer distribution in the current file
        mid = calculate_mid(ic_list, ave_left_ix, ave_right_ix)

        # store mass isotopomer distribution inside mids variable
        mids.set_values(mid, file_name)
            
    else:

        # print 'file:', file_name, 'does not have at least one mass isotopomer intensity greater than', noise 
        mid = [0] * mid_size

        # store mass isotopomer distribution inside mids variable
        mids.set_values(mid, file_name)

    return mids
