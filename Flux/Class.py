"""
Class to model MID data
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


class MIDS(object):
    
    """
    @summary: Object for storing MID values 

    @author: Milica Ng
    """

    def __init__(self, name, rt, ion_list, mid_size):
        
        """
        @summary: Initialise the MIDS data object

        @param name: Name of the compound
        @type name: StringType
        @param rt: Chromatographic retention time
        @type rt: FloatType
        @param ion_list: List of compound diagnostic ions
        @type ion_list: ListType
        @param mid_size: total number of masses (n+1 for M, M+1, ..., M+n)
        @type mid_size: IntType

        @author: Milica Ng
        """
        
        self.__name = name
        self.__rt = rt
        self.__ion_list = ion_list
        self.__mid_size = mid_size
        self.__values = {}

    def get_name(self):

        """
        @summary: Return the compound name

        @return: Name of the compound
        @rtype: StringType

        @author: Milica Ng
        """

        return self.__name

    def get_rt(self):

        """
        @summary: Return chromatographic retention time

        @return: Chromatographic retention time
        @rtype: FloatType

        @author: Milica Ng
        """

        return self.__rt

    def get_ion_list(self):

        """
        @summary: Return list of diagnostic ions

        @return: List of diagnostic ions
        @rtype: ListType

        @author: Milica Ng
        """

        return self.__ion_list

    def get_mid_size(self):

        """
        @summary: Return MID size

        @return: Mass isotopomer distribution size
        @rtype: IntType

        @author: Milica Ng
        """

        return self.__mid_size

    def set_values(self, mid, ion, file_num):

        """
        @summary: Set the MID values for an ion inside a particular file

        @param mid: Mass isotoomer distribution (MID) values
        @type mid: ListType
        @param ion: Diagnostic ion
        @type ion: IntType
        @param file_num: File number
        @type file_num: StringType

        @author: Milica Ng
        """

        self.__values[ion,file_num] = mid

    def write(self, out_file):

        """
        @summary: Write MID data to a file

        @param name: Name of the output file
        @type name: StringType

        @author: Milica Ng
        """

        print ' -> Writing to file ', out_file

        # write a header (name and retention time)
        fp = open(out_file, 'a')
        fp.write('\n')
        fp.write(self.__name)
        fp.write('\n')
        fp.write(str(self.__rt))
        fp.write(' secs')
        fp.write(',')
        fp.write(str(self.__rt/float(60)))
        fp.write(' mins')
        # fp.write(',')
        # fp.write(str(self.__mid_size))
        fp.write('\n')

        # write column names
        fp.write('ion')
        fp.write(',')
        fp.write('file num')
        fp.write(',')
        for m in range(0,self.__mid_size):
            fp.write('M+')
            fp.write(str(m))
            fp.write(',')            
        fp.write('\n')

        # write mid values
        keys = self.__values.keys()
        keys.sort()
        for k in keys:
            fp.write(str(k[0]))
            fp.write(',')
            fp.write(str(k[1]))
            mid = self.__values[k]
            fp.write(',') 
            mid_sum = float(sum(mid))
            for i in range(0, len(mid)):     
                fp.write(str(mid[i]/mid_sum))
                fp.write(',')
            fp.write('\n')
        fp.close()
        print '\n'


