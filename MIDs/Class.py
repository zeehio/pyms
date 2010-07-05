"""
Class to model MID data
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


class MID_table(object):
    
    """
    @summary: Object for storing MID values 

    @author: Milica Ng
    """

    def __init__(self, compound_name, rt, ion, mdv_size):
        
        """
        @summary: Initialise the MID table data object

        @param compound_name: Name of the compound
        @type compound_name: StringType
        @param rt: Chromatographic retention time
        @type rt: FloatType
        @param ion: Diagnostic ion
        @type ion: IntType
        @param mdv_size: total number of masses (n+1 for M, M+1, ..., M+n)
        @type mdv_size: IntType

        @author: Milica Ng
        """
        
        self.__compound_name = compound_name
        self.__rt = rt
        self.__ion = ion
        self.__mdv_size = mdv_size
        self.__values = {}
        self.__warnings = []
        self.__atoms = {}
        self.__fl = {}

    def is_empty(self):

        """
        @summary: Check if MID table is empty

        @return: True or False
        @rtype: BooleanType

        @author: Milica Ng
        """

        if self.__values == {}:
            empty = True
        else:
            empty = False

        return empty

    def get_compound_name(self):

        """
        @summary: Return the compound name

        @return: Name of the compound
        @rtype: StringType

        @author: Milica Ng
        """

        return self.__compound_name


    def get_rt(self):

        """
        @summary: Return chromatographic retention time

        @return: Chromatographic retention time
        @rtype: FloatType

        @author: Milica Ng
        """

        return self.__rt

    def get_ion(self):

        """
        @summary: Return diagnostic ion

        @return: Diagnostic ion
        @rtype: IntType

        @author: Milica Ng
        """

        return self.__ion

    def get_mdv_size(self):

        """
        @summary: Return MDV size

        @return: Mass isotopomer distribution vector size
        @rtype: IntType

        @author: Milica Ng
        """

        return self.__mdv_size

    def get_values(self):

        """
        @summary: Return MID table values

        @return: values
        @rtype: DictType

        @author: Milica Ng
        """

        return self.__values

    def get_atoms(self):

        """
        @summary: Return ion atom composition

        @return: atoms
        @rtype: DictType

        @author: Milica Ng
        """

        return self.__atoms


    def set_atoms(self, element, num):

        """
        @summary: Set the number of atoms for a particular element

        @param element: Element (C, O, N, H, Si or S)
        @type element: StringType
        @param num: Number of particular atoms in the fragment
        @type num: IntType

        @author: Milica Ng
        """

        self.__atoms[element] = num


    def set_values(self, mdv, file_name):

        """
        @summary: Set the MID vector values inside a particular file

        @param mdv: Mass isotopomer distribution vector values
        @type mdv: ListType
        @param file_name: File number
        @type file_name: StringType

        @author: Milica Ng
        """

        self.__values[file_name] = mdv

    def set_fl(self, fl, file_name):

        """
        @summary: Set the MID vector values inside a particular file

        @param fl: Fractionla labelling of the ion
        @type mdv: FloatType
        @param file_name: File number
        @type file_name: StringType

        @author: Milica Ng
        """

        self.__fl[file_name] = fl

    def append_warning(self, warning):

        """
        @summary: Append warning to the warnings list

        @param warning: Warning text
        @type warning: StringType

        @author: Milica Ng
        """

        self.__warnings[len(self.__warnings):] = [warning]

    def write(self, out_file):

        """
        @summary: Write MID table data to a file

        @param out_file: Name of the output file
        @type out_file: StringType

        @author: Milica Ng
        """

        # write a header (compound name and retention time)
        fp = open(out_file, 'a')
        fp.write('\n')
        fp.write('compound,')
        fp.write(self.__compound_name)
        fp.write('\n')
        fp.write('rt =,')
        fp.write(str(self.__rt))
        fp.write(',secs')
        fp.write('\n')
        fp.write('#rt =,')
        fp.write(str(self.__rt/float(60)))
        fp.write(',mins')
        fp.write('\n')
        fp.write('ion,')
        fp.write(str(self.__ion))
        fp.write('\n')
        fp.write('mdv size =,')
        fp.write(str(self.__mdv_size))

        fp.write('\n')

        # write column names
        fp.write('#,,')
        for m in range(0,self.__mdv_size):
            fp.write('M+')
            fp.write(str(m))
            fp.write(',') 
        fp.write('FL')           
        fp.write('\n')

        # write mdv values
        keys = self.__values.keys()
        keys.sort()

        for k in keys:

            # write file name
            fp.write('file name,')
            fp.write(str(k))
            fp.write(',')

            # write mass isotopomer distribution
            mdv = self.__values[k]
            fl = self.__fl[k]
            mdv_sum = float(sum(mdv))
            if mdv_sum > 0:
                for i in range(0, len(mdv)):     
                    fp.write(str(mdv[i]/mdv_sum))
                    fp.write(',')
                fp.write(str(fl))
                fp.write('\n')
            else:
                for i in range(0, len(mdv)):     
                    fp.write(str(0.0))
                    fp.write(',')
                fp.write('\n')

        # write warnings
        warning_list = self.__warnings
        fp.write('\n')
        for warning in warning_list:
            fp.write('\n')
            fp.write(warning)
        fp.write('\n')
        
        # close the file           
        fp.close()


