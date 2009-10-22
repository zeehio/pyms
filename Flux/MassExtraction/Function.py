"""
Function for LabelledData.MID.Extract
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


from pyms.GCMS.IO.ANDI.Function import ANDI_reader
from pyms.GCMS.Function import build_intensity_matrix_i


def extract_mid(file_base, file_numbers, compound, ion, rt, mid_size, win_size):


    """
    @summary: Method for extracting mass isotopomer distribution (MID)

    @param file_base: Base for file names
    @type file_base: types.StringType
    @param file_numbers: File numbers
    @type file_numbers: types.TupleType
    @param compound: Name of the compound
    @type compound: types.StringType
    @param ion: m/z value of the fragement ion
    @type ion: types.IntType
    @param rt: Retention time of the compound
    @type rt: types.FloatType
    @param mid_size: Size of the mass distribution vector
    @type mid_size: types.IntType
    @param win_size: Size of the window for local maximum search
    @type win_size: types.FloatType

    @return: None
    @rtype: types.NoneType
    
    @author: Milica Ng
    """

    print '\n'
    print 'Compound:', compound

    for file_num in file_numbers:

        left_boundary_sum = 0
        right_boundary_sum = 0
        sum_count = 0
        andi_file = file_base+str(file_num)+".CDF"

        # load raw data
        data = ANDI_reader(andi_file)

        # create intensity matrix
        im = build_intensity_matrix_i(data)

        for mz in range(ion, ion+mid_size):
            # get ion chromatogram at specified m/z
            ic = im.get_ic_at_mass(mz)

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
            peak_apex_rt = ic.get_time_at_index(peak_apex_ix)

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
            left_boundary_rt = ic.get_time_at_index(left_boundary_ix)

            # find right_boundary (first local minimum to the right of peak apex)
            right_boundary_int = peak_apex_int
            ix = peak_apex_ix
            while ix <= len(ic):
                if ic.get_intensity_at_index(ix) > right_boundary_int:
                    break
                else:
                    right_boundary_int = ic.get_intensity_at_index(ix)
                    ix = ix + 1
            right_boundary_ix = ix - 1
            right_boundary_rt = ic.get_time_at_index(right_boundary_ix)

            # print left_boundary_ix, peak_apex_ix, right_boundary_ix, peak_apex_int

            # find average of left and right boundary index for the file being processed
            if peak_apex_int > 4000: # todo: estimate 4000 from the file
                left_boundary_sum = left_boundary_sum + left_boundary_ix
                right_boundary_sum = right_boundary_sum + right_boundary_ix
                sum_count = sum_count + 1
        if sum_count > 0:
            ave_left_ix = left_boundary_sum/sum_count
            ave_right_ix = right_boundary_sum/sum_count
            # print ' ave_left_ix' , ave_left_ix
            # print ' ave_right_ix', ave_right_ix
            # print ' sum_count', sum_count
            print 'Mass Isotopomer Distribution (MID):',
            for mz in range(ion, ion+mid_size):
                # get ion chromatogram at specified m/z
                ic = im.get_ic_at_mass(mz)
                area = 0
                # integrate (sum intensity from average left_boundary_ix to average right_boundary_ix 
                for ix in range(ave_left_ix, ave_right_ix+1):
                    area = area + ic.get_intensity_at_index(ix)
                print area,
            print '\n'
        else:
            print 'file:', andi_file, 'does not have at least one mass isotopomer intensity greater than 4000'
