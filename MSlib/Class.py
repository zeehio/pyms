"""
Class to create mass spectral library
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
sys.path.append("/x/PyMS/")

class msrecord(object):
    
    """
    @summary: Models a mass spectral record

    @author: Saravanan Dayalan
    """

    def __init__(self):
        
        """     
        @summary: Initialises an instance of a mass spectral record
        """
        
        self.name = ''
        self.regno = ''
        self.mi = {}


    def set(self, name, regno, mi):
        
        """
        @summary: Assigns name, CAS registration number
        and mass spectrum to a mass spectral record
        
        @param name: Name of the compound 
        @type name: StringType
        
        @param regno: CAS registration number of the compound
        @type regno: StringType

        @param mi: Mass spectral values (mass, intensity) of
        the compound
        @type mi: listType
        """
        
        self.name = name
        self.regno = regno
        self.mi = mi


class mslib(object):
    
    """
    @summary: Models a mass spectral library

    @author: Saravanan Dayalan
    """

    def __init__(self):
        
        """     
        @summary: Initialises an instance of a mass spectral library
        """
                
        self.msrecord_list = []

    def addrecord(self, msrecord):
        
        """
        @summary: Adds a mass spectral record to the
        mass spectral library       
        
        @param msrecord: An instance of a mass spectral record
        @type msrecord: pyms.MSlib.Class.msrecord
        """
                
        self.msrecord_list.append(msrecord)

