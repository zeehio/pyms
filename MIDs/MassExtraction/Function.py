"""
Functions for Flux.MassExtract
"""


 #############################################################################
 #                                                                           #
 #    PyMS software for processing of metabolomic mass-spectrometry data     #
 #    Copyright (C) 2005-2010 Vladimir Likic                                 #
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


from pyms.MIDs.Class import MID
from pyms.Noise.SavitzkyGolay import savitzky_golay
from pyms.Baseline.TopHat import tophat

# for plotting only, delete later!
from pylab import * # used by plotfile, savefig

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
        # print 'ic.get_intensity_at_index(ix)', ic.get_intensity_at_index(ix), 'left_boundary_int', left_boundary_int
        
        if ic.get_intensity_at_index(ix) > left_boundary_int or 0 == left_boundary_int:
            # print 'broke on left'
            # print 'ic.get_intensity_at_index(ix+1)', ic.get_intensity_at_index(ix+1),'left_boundary_int', left_boundary_int
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
        # print 'ic.get_intensity_at_index(ix)', ic.get_intensity_at_index(ix),'right_boundary_int', right_boundary_int
        if ic.get_intensity_at_index(ix) > right_boundary_int or 0 == right_boundary_int:
            # print ' broke on right'
            # print 'ic.get_intensity_at_index(ix+1)', ic.get_intensity_at_index(ix+1), 'right_boundary_int', right_boundary_int
            break
        else:
            right_boundary_int = ic.get_intensity_at_index(ix)
            ix = ix + 1
    right_boundary_ix = ix - 1

    return right_boundary_ix


def calculate_mdv(ic_list, ave_left_ix, ave_right_ix):

    """
    @summary: Calculate mass isotopomer distribution vector

    @param ic_list: List of IonChromatograms
    @type ic_list: ListType
    @param ave_left_ix: Average left boundary index for the ions
    @type ave_left_ix: IntType
    @param ave_right_ix: Average right boundary index for the ions
    @type ave_right_ix: IntType
    
    @return: Mass isotopomer distribution vector
    @rtype: ListType
    
    @author: Milica Ng
    """

    # calculate mass isotopomer distribution
    mdv = []
    for ic in ic_list:
        # get ion chromatogram at specified m/z
        # ic = im.get_ic_at_mass(mz)
        area = 0
        # integrate (sum intensity from average left_boundary_ix to average right_boundary_ix 
        for ix in range(ave_left_ix, ave_right_ix+1):
            area = area + ic.get_intensity_at_index(ix)
        mdv.append(area)

    return mdv

def plot_ics(ic_list, ave_left_ix, ave_right_ix, file_name, compound, noise):

    """
    @summary: For testing purposes only! To be deleted.
    
    @author: Milica Ng
    """

    count = 0
    for ic in ic_list:
        col1 = []
        col2 = []

        for ix in range(ave_left_ix-11,ave_right_ix+11):

            if (0 < ix) and (ix < len(ic)):
                col1.append(ic.get_time_at_index(ix))
                col2.append(ic.get_intensity_at_index(ix))

        plot_file = "plots/"+compound+"-"+str(file_name)+"-"+str(count)+".png"
        count = count+1        
        plot(col1,col2) 
        vlines(ic.get_time_at_index(ave_left_ix), 0, noise)
        vlines(ic.get_time_at_index(ave_right_ix), 0, noise)
        title(plot_file)
        savefig(plot_file)
        close()


def extract_mid(mids, file_name, im, win_size, noise):

    """
    @summary: Method for extracting mass isotopomer distribution (MID)

    @param file_name: File number
    @type file_name: StringType
    @param im: Intensity matrix
    @type im: pyms.GCMS.Class.IntensityMatrix
    @param mids: Mass isotopomer distribution data
    @type mids: pyms.Flux.Class.MIDS
    @param mdv_size: Size of the mass distribution vector
    @type mdv_size: IntType
    @param win_size: Size of the window for local maximum search
    @type win_size: FloatType
    @param noise: Noise threshold below which signal is assumed unreliable
    @type win_size: IntType
    
    @author: Milica Ng
    """
    
    # get compound name, rt, ion and mdv size
    compound = mids.get_compound_name()
    rt = mids.get_rt()
    ion = mids.get_ion()
    mdv_size = mids.get_mdv_size()

    # initialise loop variables
    left_boundary_sum = 0
    right_boundary_sum = 0
    sum_count = 0
    ic_list = []
    peak_apex_time_list = []
    peak_apex_intensity_list = []

    # loop over required number of mass isotopomers
    for mz in range(ion, ion+mdv_size):

        # get ion chromatogram at current m/z
        # todo: check if mass is out of range, raise an error!
        ic = im.get_ic_at_mass(mz) 

        # smooth noise + correct baseline
        ic_smooth = savitzky_golay(ic)
        ic_bc = tophat(ic_smooth, struct="1.5m")
        ic = ic_bc

        # store ic's for calculate_mdv later
        ic_list.append(ic)

        # get peak apex index and intensity
        peak_apex_ix = peak_apex(ic, rt, win_size)
        peak_apex_int = ic.get_intensity_at_index(peak_apex_ix)

        # store peak apex retention time & intensity for later checks
        peak_apex_time_list.append(ic.get_time_at_index(peak_apex_ix))
        peak_apex_intensity_list.append(ic.get_intensity_at_index(peak_apex_ix))

        # get indices for left/right boundary
        left_boundary_ix = left_boundary(ic, peak_apex_ix, peak_apex_int)
        right_boundary_ix = right_boundary(ic, peak_apex_ix, peak_apex_int)

        if peak_apex_int > noise: 

            # raise warning if current peak apex retention time is within half a window of previous ones
            for time in peak_apex_time_list:
                if (abs(ic.get_time_at_index(peak_apex_ix) - time)) > (win_size/2):
                     warning = 'Warning: '+compound+' '+str(file_name)+' '+str(mz)+' Peak apex rt shift > '+str(win_size/2)+'secs'
                     mids.append_warning(warning)
                     break
            # raise warning if peak apex is smaller than left or right intensity
            if ic.get_intensity_at_index(peak_apex_ix) < ic.get_intensity_at_index(peak_apex_ix-1):
                warning = 'Warning: '+compound+' '+str(file_name)+' '+str(mz)+' An intensity > peak apex on the LEFT'
                mids.append_warning(warning)
            elif ic.get_intensity_at_index(peak_apex_ix) < ic.get_intensity_at_index(peak_apex_ix+1):
                warning = 'Warning: '+compound+' '+str(file_name)+' '+str(mz)+' An intensity > peak apex on the RIGHT'
                mids.append_warning(warning)
            # raise warning if peak apex equals both neighbouring intensities 
            elif (ic.get_intensity_at_index(peak_apex_ix) == ic.get_intensity_at_index(peak_apex_ix-1)) and (ic.get_intensity_at_index(peak_apex_ix) == ic.get_intensity_at_index(peak_apex_ix+1)):
                warning = 'Warning: '+compound+' '+str(file_name)+' '+str(mz)+' Peak apex intensity EQUALS neighbouring left and right intensity'
                mids.append_warning(warning)

            # sum left and right boundary indices above noise in the current file
            left_boundary_sum = left_boundary_sum + left_boundary_ix
            right_boundary_sum = right_boundary_sum + right_boundary_ix
            sum_count = sum_count + 1

            # collect left and right boundaries for standard deviation later
            # left_boundary_ix_list.append(left_boundary_ix)
            # right_boundary_ix_list.append(right_boundary_ix)

    if sum_count > 0:


        # raise warning if peak boundaries belonging to a single peak were used as average
        if sum_count == 1:
            warning = 'Warning: '+compound+' '+str(file_name)+' Only one ion peak was larger than '+str(noise)+'!'
            mids.append_warning(warning)

        # find average of left and right boundary index in the current file
        ave_left_ix = left_boundary_sum/sum_count
        ave_right_ix = right_boundary_sum/sum_count

        # raise warning if rt is not between ave left and right boundary
        if (rt < ave_left_ix) and (rt > ave_right_ix):
            warning = 'Warning: '+compound+' '+str(file_name)+' Retention time is outside the integration interval'
            mids.append_warning(warning)

        # raise warning if an intensity larger than noise is found where peak apex is smaller than noise
        j = 0
        for peak_apex_intensity in peak_apex_intensity_list:
             if peak_apex_intensity < noise:
                 for ix in range(ave_left_ix, ave_right_ix):
                     if ic_list[j].get_intensity_at_index(ix) > noise:
                         warning = 'Warning: '+compound+' '+str(file_name)+' '+str(ion+j)+' Peak apex detected  as lower than '+str(noise)+' however an intensity integrated was larger than '+str(noise)
                         mids.append_warning(warning)
                         break
             j = j+1

        # plot for testing purposes only. to be deleted! 
        # plot_ics(ic_list, ave_left_ix, ave_right_ix, file_name, compound, noise)

        # calculate mass isotopomer distribution in the current file
        mdv = calculate_mdv(ic_list, ave_left_ix, ave_right_ix)
            
    else:
        mdv = [0] * mdv_size

    # store mass isotopomer distribution inside mids variable
    mids.set_values(mdv, file_name)


