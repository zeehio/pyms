"""
Models a GC-MS experiment represented by a list of signal peaks
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
from pyms.Utils.Utils import is_list, is_str
from pyms import Peak

class Experiment:

    """
    @summary: Models an experiment object

    @author: Vladimir Likic
    @author: Andrew Isaac
    """

    def __init__(self, name, peak_list):

        """
        @summary: Models an experiment

        @param name: Unique identifier for the experiment
        @type name: StringType
        @param peak_list: A list of peak objects
        @type peak_list: ListType
        """

        if not is_str(name):
            error("'name' must be a string")
        if not is_list(peak_list):
            error("'peak_list' must be a list")
        if not len(peak_list) > 0 and not isinstance(peak_list[0], Peak):
            error("'peak_list' must be a list of Peak objects")

        self.__name = name
        self.__peak_list = peak_list

    def get_name(self):

        """
        @summary: Returns the name of the experiment

        @return: The name of the experiment
        @rtype:  StringType
        """

        return self.__name

    def get_peak_list(self):

        """
        @summary: Returns the peak list

        @return: A list of peak objects
        @rtype: ListType
        """

        return self.__peak_list
