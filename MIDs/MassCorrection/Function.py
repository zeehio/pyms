"""
Functions for MIDs.MassCorrect
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

import numpy

from pyms.Utils.Error import error
from pyms.Utils.Utils import is_int, is_positive_int, is_list, \
    is_list_of_dec_nums, is_array, is_number
from Constants import nsi


def correct_mdv(mdv, atoms, f_unlabelled):

    """
    @summary: Correct MID for natural isotopic abundance

    @param mdv: Mass isotopomer distribution vector
    @type mdv: types.ListType
    @param atoms: Fragment's elemental composition
    @type atoms: types.DictType
    @param f_unlabelled: Fractional labelling.
    @type f_unlabelled: types.FloatType
    
    @author: Milica Ng
    """
    # calculate n, this can work depending on the content of the  input file
    n = len(mdv)-1

    # the overall correction matrix
    c_corr = overall_correction_matrix(n, mdv, atoms)
    #print('\n Overall correction matrix:')
    #print(c_corr)

    # the exclusive mass isotope distribution of the carbon skeleton, mdv_alpha_star
    mdv_alpha_star = c_mass_isotope_distr(mdv, c_corr)
    #print('\n Normalised mdv alpha star:')
    #print(mdv_alpha_star)

    # correction for unlabelled biomass
    mdv_aa = corr_unlabelled(n, mdv_alpha_star, f_unlabelled)
    #print('\n mdv_aa:')
    #print(mdv_aa)

    # For [U-13C]glucose experiments, the fractional labelling (FL) of the different
    # fragments should be equal to the labelling content of the input substrate 
    # (i.e. when 20% [U-13C]glucose is used, the FL should be 0.2 for all fragments).
    fl = fract_labelling(n, mdv_aa)
    #print('\n Fractional labelling FL: %s' % ( fl ))

    # convert to a list of lists 
    mdv_aa = mdv_aa.tolist()
    # reduce to a list
    mdv_aa = reduce(list.__add__, mdv_aa)

    return (mdv_aa, fl)

def fract_labelling(n, mdv_aa):

    """
    @summary: Calculates the fractional labelling of the fragment.

    @param n: Number of fragment's C atoms which may contain exogenous (non
        natural abundance) 13C (e.g. only amino acid ones)
    @type n: types.IntType
    @param mdv_aa: Mass distribution vector.
    @type mdv_aa: numpy.ndarray

    @return: Fractional labelling.
    @rtype: types.FloatType
    
    @author: Milica Ng
    """

    # check the arguments
    if not is_positive_int(n):
        error("'n' must be an integer greater than zero")
    if not is_array(mdv_aa):
        error("'mdv_aa' be numpy.ndarray")

    mdv_aa = numpy.matrix(mdv_aa)
    fl = (range(0,(n+1)) * mdv_aa) / (n * numpy.sum(mdv_aa))
    return float(fl)

def corr_unlabelled(n, mdv_alpha_star, f_unlabelled):

    """
    @summary: Corrects the mass distribution vector for fractional labelling.

    @param n: Number of fragment's C atoms which may contain exogenous (non
        natural abundance) 13C (e.g. only amino acid ones)
    @type n: types.IntType
    @param mdv_alpha_star: Mass distribution vector.
    @type mdv_alpha_star: numpy.ndarray
    @param f_unlabelled: Percentage unlabelled biomass.
    @type f_unlabelled: types.IntType or types.FloatType

    @return: Mass distribution vector.
    @rtype: numpy.ndarray

    @author: Milica Ng
    """

    # check the arguments
    if not is_positive_int(n):
        error("'n' must be an integer greater than zero")
    if not is_array(mdv_alpha_star):
        error("'mdv_alpha_star' must be numpy.ndarray")
    if not is_number(f_unlabelled):
        error("'f_unlabelled' must be a number")

    mdv_unlabelled = mass_dist_vector(n, n, nsi['c'])
    mdv_aa = (mdv_alpha_star - f_unlabelled * mdv_unlabelled)/(1 - f_unlabelled)

    return mdv_aa

def overall_correction_matrix(n, mdv, atoms):

    """
    @summary: Calculates the overall correction matrix.

    @param n: Number of fragment's C atoms which may contain exogenous (non
        natural abundance) 13C (e.g. only amino acid ones)
    @type n: types.IntType
    @param mdv: Mass distribution vector.
    @type mdv: types.ListType
    @param atoms: Number of C, O, N, H, Si, and S atoms in the fragment
        (note: the number of C atoms excludes C atoms which may 
        contain exogenous (non natural abundance) 13C)
    @type atoms: types.DictType
    @param nsi: Contains abundance values of natural stable isotopes.
    @type nsi: types.ListType

    @return: Overall correction matrix.
    @rtype: numpy.ndarray
    
    @author: Milica Ng
    """

    # check the arguments
    if not is_positive_int(n):
        error("'n' must be an integer greater than zero")
    if not (is_list(mdv) or not is_list(nsi)):
        error("'mdv' and 'nsi' must be types.ListType")

    atom_symbols = atoms.keys()
    c_corr = numpy.eye(n+1) 
    for a in atom_symbols:
        # print '\n Calculating %s correction matrix' % ( a )
        m_corr = correction_matrix(n, atoms[a], nsi[a])
        c_corr = numpy.dot(c_corr,m_corr)
    # print '\n Calculated overall correction matrix.'
    return c_corr

def correction_matrix(n, num_a, nsil):

    """
    @summary: Calculates a correction matrix.

    @param n: Number of fragment's C atoms which may contain exogenous (non
        natural abundance) 13C (e.g. only amino acid ones)
    @type n: types.IntType
    @param num_a: Number of C, O, N, H, Si, or S atoms in the fragment (note:
        the number of C atoms excludes C atoms which may contain exogenous
        (non natural abundance) 13C).
    @type num_a: types.IntType
    @param nsil: Contains abundance values of natural stable isotopes.
    @type nsil: types.ListType

    @return: Correction matrix.
    @rtype: numpy.ndarray
    
    @author: Milica Ng
    """

    # check the arguments
    if not is_positive_int(n):
        error("'n' must be an integer greater than zero")
    if not is_int(num_a):
        error("'num_a' must be an integer")
    if not is_list_of_dec_nums(nsil):
        error("'nsil' must be a non-empty list of decimal numbers")

    if not num_a == 0:
        c_corr_a = numpy.zeros((n+1,n+1))
        # calculate the first column of the correction matrix
        tmp = mass_dist_vector(n, num_a, nsil)
        c_corr_a[:,0]  = numpy.reshape(tmp,(1,len(tmp)))
        # calculate the rest of the matrix from the first column
        for i in range(1,(n+1)):
            for j in range(1,(n+1)):
                c_corr_a[i,j] = c_corr_a[i-1,j-1]
    else:
        c_corr_a = numpy.identity(n+1)
    return c_corr_a

def c_mass_isotope_distr(mdv, c_corr):

    """
    @summary: Calculates the mass isotope distribution of the carbon skeleton.

    @param mdv: Mass distribution vector.
    @type mdv: types.ListType
    @param c_corr: The overall correction matrix.
    @type c_corr: numpy.ndarray
    
    @return: The exclusive mass isotope distribution of the carbon skeleton.
    @rtype: numpy.ndarray
    
    @author: Milica Ng
    """

    # check the arguments
    if not is_list(mdv):
        error("'mdv' must be types.ListType")
    if not is_array(c_corr):
        error("'c_corr' must be numpy.ndarray")

    mdv = numpy.array(mdv, float)
    mdv_alpha = mdv/sum(mdv)
    mdv_alpha = numpy.matrix(mdv_alpha)
    mdv_alpha = mdv_alpha.T
    c_corr = numpy.matrix(c_corr)
    mdv_alpha_star = (c_corr.I) * mdv_alpha
    mdv_alpha_star = numpy.array(mdv_alpha_star, float)
    mdv_alpha_star = mdv_alpha_star/sum(mdv_alpha_star)
    return mdv_alpha_star

def mass_dist_vector (n, num_a, nsil):

    """
    @summary: Calculates mass distribution vector.

    @param n: Number of fragment's C atoms which may contain exogenous (non
        natural abundance) 13C (e.g. only amino acid ones)
    @type n: types.IntType
    @param num_a: Number of C, O, N, H, Si, or S atoms in the fragment
        (note: the number of C atoms excludes C atoms which may contain
        exogenous (non natural abundance) 13C).
    @type num_a: types.IntType
    @param nsil: Contains abundance values of natural stable isotopes.
    @type nsil: types.ListType of types.FloatType

    @return: Mass distribution vector.
    @rtype: numpy.ndarray

    @author: Milica Ng
    """

    # check the arguments
    if not is_positive_int(n):
        error("'n' must be an integer greater than zero")
    if not is_positive_int(num_a):
        error("'n' must be an integer greater than zero")
    if not is_list_of_dec_nums(nsil):
        error("'nsil' must be a non-empty list of decimal numbers")

    mdv = numpy.zeros(((n + 1),1))
    #   calculate all mass combinations where value 0 corresponds to m0,
    #     1 corresponds to m1, 2 corresponds to m2, ..., N corresponds to mN
    m_combs = []
    m_combs = combinations_with_repetition(numpy.size(range(0,(len(nsil)))), num_a) - 1
    # sort m_combs in order of row sums
    m_combs = m_combs.tolist()
    m_combs.sort(lambda x, y: cmp(sum(x),sum(y)))
    m_combs = numpy.array(m_combs)
    m = 0  # initial mass value
    count = 0  # current row in m_comb
    for i in range(0,n+1):
        # loop until the value of m changes (e.g m0 -> m1)
        while m == sum(m_combs[count,]):
            # calculate v vector
            v = numpy.zeros((len(nsil),1))
            for j in range(0,num_a):
                v[(m_combs[count,j])] = v[(m_combs[count,j])] + 1
            # calculate isotopolog abundance (part_1 * part_2) and add any
            # previous value if more than one combination corresponds to mN
            part_1 = fact(int(sum(v)))
            part_2 = 1
            for k in range(0,(len(nsil))):
                part_2 = part_2 * nsil[k]**v[k]/fact(int(v[k]))
            mdv[i] = mdv[i] + part_1 * part_2
            count = count + 1
            if (count > (numpy.size(m_combs,0)) - 1):  # check if finished
                break                              # (some entries might 0)
        if (((i + 1) > (numpy.size(m_combs,0)) - 1) or (count > (numpy.size(m_combs,0)) - 1)):
            break                              # (some entries might be 0)
        m = m + 1
    return mdv

def fact(integer):
    # from function import fact
    factorial = 1
    for factor in range(2, integer+1):
        factorial *= factor
    return factorial

def combinations_with_repetition(n, k):

    """
    @summary: Calculates indices for picking k elements (with replacement,
        order matters) from the set with n elements.

    @param n: Number of elements in the set
    @type n: types.IntType
    @param k: Number of picks
    @type k: types.IntType

    @return: Combinations with repetition indices (base 1).
    @rtype: numpy.ndarray

    @author: Milica Ng
    """

    # check the arguments
    if not is_positive_int(n):
        error("'n' must be an integer greater than zero")
    if not is_positive_int(k):
        error("'k' must be an integer greater than zero")

    if (k==1):
        indices = numpy.arange(1,n+1)[:, numpy.newaxis]
        return indices
    if (n==1):
        indices = numpy.ones((1,k,))
        return indices
    indices = []
    for z in range(1,n+1):
        next_indices = combinations_with_repetition(n+1-z,k-1)
        if indices == []:
            indices = numpy.hstack((z*numpy.ones((numpy.size(next_indices,0),1)), next_indices+z-1))
        else:
            tmp = numpy.hstack((z*numpy.ones((numpy.size(next_indices,0),1)), next_indices+z-1))
            indices = numpy.vstack((indices, tmp))
    return indices
