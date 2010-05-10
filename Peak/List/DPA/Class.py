"""
Classes for peak alignment by dynamic programming
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

import copy

import numpy, Pycluster

from pyms.Utils.Error import error, stop
from pyms.Experiment.Class import Experiment
from pyms.GCMS.Class import MassSpectrum
from pyms.Peak.Class import Peak
from pyms.Peak.List.Function import composite_peak

import Function

# If psyco is installed, use it to speed up running time
try:
    import psyco
    psyco.full()
except:
    pass

class Alignment(object):

    """
    @summary: Models an alignment of peak lists

    @author: Woon Wai Keen
    @author: Vladimir Likic
    """

    def __init__(self, expr):

        """
        @param expr: The experiment to be converted into an alignment object
        @type expr: pyms.Experiment.Class.Experiment

        @author: Woon Wai Keen
        @author: Qiao Wang
        @author: Vladimir Likic
        """
        if expr == None:
            self.peakpos = []
            self.peakalgt = []
            self.expr_code =  []
            self.similarity = None
        else:
            if not isinstance(expr, Experiment):
                error("'expr' must be an Experiment object")
            for peak in expr.get_peak_list():
                if peak.get_area() == None or peak.get_area() <= 0:
                    error("All peaks must have an area for alignment")
            self.peakpos = [ copy.deepcopy(expr.get_peak_list()) ]
            self.peakalgt = numpy.transpose(self.peakpos)
            self.expr_code =  [ expr.get_expr_code() ]
            self.similarity = None

    def __len__(self):

        """
        @summary: Returns the length of the alignment, defined as the number of
        peak positions in the alignment

        @author: Qiao Wang
        @author: Vladimir Likic
        """

        return len(self.peakalgt)

    def filter_min_peaks(self, min_peaks):

        """
        @summary: Filters alignment positions that have less peaks than 'min_peaks'

        This function is useful only for within state alignment.

        @param min_peaks: Minimum number of peaks required for the alignment
            position to survive filtering
        @type min_peaks: IntType

        @author: Qiao Wang
        """

        filtered_list=[]

        for pos in range(len(self.peakalgt)):
            if len(filter(None,self.peakalgt[pos])) >= min_peaks:
                filtered_list.append(self.peakalgt[pos])

        self.peakalgt = filtered_list
        self.peakpos = numpy.transpose(self.peakalgt)

    def write_csv(self, rt_file_name, area_file_name, minutes=True):

        """
        @summary: Writes the alignment to CSV files

        This function writes two files: one containing the alignment of peak
        retention times and the other containing the alignment of peak areas.

        @param rt_file_name: The name for the retention time alignment file
        @type rt_file_name: StringType
        @param area_file_name: The name for the areas alignment file
        @type area_file_name: StringType
        @param minutes: An optional indicator whether to save retention times
            in minutes. If False, retention time will be saved in seconds
        @type minutes: BooleanType

        @author: Woon Wai Keen
        @author: Andrew Isaac
        @author: Vladimir Likic
        """

        try:
            fp1 = open(rt_file_name, "w")
            fp2 = open(area_file_name, "w")
        except IOError:
            error("Cannot open output file for writing")

        # create header
        header = '"UID","RTavg"'
        for item in self.expr_code:
            expr_code = ( '"%s"' % item )
            header = header + "," + expr_code
        header = header + "\n"

        # write headers
        fp1.write(header)
        fp2.write(header)

        # for each alignment position write alignment's peak and area
        for peak_idx in range(len(self.peakpos[0])):

            rts = []
            areas = []
            new_peak_list = []
            avgrt = 0
            countrt = 0

            for align_idx in range(len(self.peakpos)):

                peak = self.peakpos[align_idx][peak_idx]

                if peak is not None:

                    if minutes:
                        rt = peak.get_rt()/60.0
                    else:
                        rt = peak.get_rt()

                    rts.append(rt)
                    areas.append(peak.get_area())
                    new_peak_list.append(peak)

                    avgrt = avgrt + rt
                    countrt = countrt + 1
                else:
                    rts.append(None)
                    areas.append(None)

            if countrt > 0:
                avgrt = avgrt/countrt

            compo_peak = composite_peak(new_peak_list, minutes)
            peak_UID = compo_peak.get_UID()
            peak_UID_string = ( '"%s"' % peak_UID)

            # write to retention times file
            fp1.write(peak_UID_string)
            fp1.write(",%.3f" % avgrt)
            for rt in rts:
                if rt == None:
                    fp1.write(",NA")
                else:
                    fp1.write(",%.3f" % rt)
            fp1.write("\n")

            # write to peak areas file
            fp2.write(peak_UID_string)
            fp2.write(",%.3f" % avgrt)
            for area in areas:
                if area == None:
                    fp2.write(",NA")
                else:
                    fp2.write(",%.4f" % area)
            fp2.write("\n")

        fp1.close()
        fp2.close()

    def aligned_peaks(self, minutes=False):

        """
        @summary: Returns a list of Peak objects where each peak
            has the combined spectra and average retention time
            of all peaks that aligned.

        @param minutes: An optional indicator of whether retention
            times are in minutes. If False, retention time are in
            seconds
        @type minutes: BooleanType

        @return: A list of composite peaks based on the alignment.
        @rtype: ListType

        @author: Andrew Isaac
        """

        # for all peaks found
        peak_list = []
        for peak_idx in range(len(self.peakpos[0])):
            # get aligned peaks, ignore missing
            new_peak_list = []
            for align_idx in range(len(self.peakpos)):
                peak = self.peakpos[align_idx][peak_idx]
                if peak is not None:
                    new_peak_list.append(peak)
            #create composite
            new_peak = composite_peak(new_peak_list, minutes)
            peak_list.append(new_peak)

        return peak_list

class PairwiseAlignment(object):

    """
    @summary: Models pairwise alignment of alignments

    @author: Woon Wai Keen
    @author: Vladimir Likic
    """

    def __init__(self, algts, D, gap):

        """
        @param algts: A list of alignments
        @type algts: ListType
        @param D: Retention time tolerance parameter for pairwise alignments
        @type D: FloatType
        @param gap: Gap parameter for pairwise alignments
        @type gap: FloatType

        @author: Woon Wai Keen
        @author: Vladimir Likic
        """
        self.algts = algts
        self.D = D
        self.gap = gap

        self.sim_matrix = self._sim_matrix(algts, D, gap)
        self.dist_matrix = self._dist_matrix(self.sim_matrix)
        self.tree = self._guide_tree(self.dist_matrix)

    def _sim_matrix(self, algts, D, gap):

        """
        @summary: Calculates the similarity matrix for the set of alignments

        @param algts: A list of alignments
        @type algts: ListType
        @param D: Retention time tolerance parameter for pairwise alignments
        @type D: FloatType
        @param gap: Gap parameter for pairwise alignments
        @type gap: FloatType

        @author: Woon Wai Keen
        @author: Vladimir Likic
        """

        n = len(algts)

        total_n = n * (n - 1) / 2

        print " Calculating pairwise alignments for %d alignments (D=%.2f, gap=%.2f)" % \
                (n, D, gap)

        sim_matrix = numpy.zeros((n,n), dtype='f')

        for i in range(n - 1):
            for j in range(i + 1, n):
                ma = Function.align(algts[i], algts[j], D, gap)
                sim_matrix[i,j] = sim_matrix[j,i] = ma.similarity
                total_n = total_n - 1
                print " -> %d pairs remaining" % total_n

        return sim_matrix

    def _dist_matrix(self, sim_matrix):

        """
        @summary: Converts similarity matrix into a distance matrix

        @param sim_matrix: The similarity matrix
        @type sim_matrix: numpy.ndarray

        @author: Woon Wai Keen
        @author: Vladimir Likic
        """

        # change similarity matrix entries (i,j) to max{matrix}-(i,j)
        sim_max = numpy.max(numpy.ravel(sim_matrix))
        dist_matrix = sim_max - sim_matrix

        # set diagonal elements of the similarity matrix to zero
        for i in range(len(dist_matrix)):
            dist_matrix[i,i] = 0

        return dist_matrix

    def _guide_tree(self, dist_matrix):

        """
        @summary: Build a guide tree from the distance matrix

        @param dist_matrix: The distance matrix
        @type dist_matrix: numpy.ndarray
        @return: Pycluster similarity tree
        @rtype: Pycluster.cluster.Tree

        @author: Woon Wai Keen
        @author: Vladimir Likic
        """

        n = len(dist_matrix)

        print " -> Clustering %d pairwise alignments." % (n*(n-1)),
        tree = Pycluster.treecluster(distancematrix=dist_matrix, method='a')
        print "Done"

        return tree

