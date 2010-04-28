"""
Provides functions for simulation of GCMS data
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
import math
import sys

sys.path.append("/x/PyMS/")
from pyms.GCMS.Class import IonChromatogram


def gcms_sim(n_scan, n_mz, period, t1, peak_list):
    
    """
    @summary: Simulator of GCMS data
    
    @param n_scan: the number of scans
    @type n_scan: intType
    
    @param n_mz: the number of m/z channels
    @type n_mz: intType
    
    @param period: the sampling interval
    @type period: floatType
    
    @param t1: The start time for sampling
    @type t1: floatType
    
    @param peak_list: A list of peaks
    @type peak_list: list of pyms.Peak.Class.Peak
    
    @author: Sean O'Callaghan
    """
    
    # initialise a 2D numpy array for intensity matrix
    im = numpy.zeros((n_scan, n_mz), 'd') 
    
    
    for peak in peak_list:
        print "-", 
        index = int((peak.get_rt() - t1)/period) 
        height = sum(peak.get_mass_spectrum().mass_spec)
        # standard deviation = area/(height * sqrt(2/pi)) 
        sigma = peak.get_area() / (height * (math.sqrt(2*math.pi)))
        print "width", sigma
        for i in range(len(peak.ms.mass_list)):
            ion_height = peak.ms.mass_spec[i]
            ic = chromatogram(n_scan, index, sigma, ion_height)
            im[:, i] += ic
            
    return im
        



def chromatogram(n_scan, x_zero, sigma, peak_scale):

    """
    @summary: Returns a simulated ion chromatogram of a pure component
              The ion chromatogram contains a single gaussian peak.
    
    @param n_scan: the number of scans
    @type n_scan: intType
    
    @param x_zero: The apex of the peak
    @type x_zero: intType
    
    @param sigma: The standard deviation of the distribution
    @type sigma: floatType
    
    @param: peak_scale: the intensity of the peak at the apex
    @type peak_scale: floatType
    
    @author: Sean O'Callaghan
    """

   
    ic = numpy.zeros((n_scan), 'd')                     
    
    for i in range(n_scan):
	x = float(i)
	
        ic[i] = gaussian(x,x_zero,sigma,peak_scale)

    return ic

def gaussian(point, mean, sigma, scale):
    """
    @summary: calculates a point on a gaussian density function
    
    f = s*exp(-((x-x0)^2)/(2*w^2));
    
    @param point: The point currently being computed
    @type point: floatType
    
    @param mean: The apex of the peak
    @type mean: intType
    
    @param width: The standard deviation of the gaussian
    @type width: floatType
     
    @param scale: The height of the apex
    @type scale: floatType
    
    @author: Sean O'Callaghan
    """
     
    return scale*math.exp((-(point-mean)**2)/(2*(sigma**2)))


def add_const_noise(ic, scale):
    """
    @summary: adds noise drawn from a normal distribution 
    with constant scale to an ion chromatogram
    
    @param ic: The ion Chromatogram
    @type ic: pyms.GCMS.IonChromatogram
    
    @param scale: The scale of the normal distribution
    @type scale: intType
    
    @author: Sean O'Callaghan
    """
    
    noise = numpy.random.normal(0.0, scale, (len(ic)))
     
    ic_array_with_noise = ic.get_intensity_array() + noise
    time_list = ic.get_time_list()
    
    ic_with_noise = IonChromatogram(ic_array_with_noise, time_list)
    
    return ic_with_noise
