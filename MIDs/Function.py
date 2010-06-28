"""
Provides helper functions for MID processing
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

from pyms.Utils.Error import error
from pyms.MIDs.Class import MID_table
from pyms.Utils.IO import file_lines
from pyms.Utils.Time import time_str_secs

def parse_ion_defs(in_file):

    """
    @summary: Read ion definitions and return as a list of MID_table objects

    @param in_file: The name of the file containing ion definitions
    @type in_file: StringType

    @return: Empty MID_table object list with compound names, 
        retention times, diagnostic ions and mdv sizes set
    @rtype: ListType

    @author: Milica Ng
    @author: Vladimir Likic
    """

    lines = file_lines(in_file, filter=True)
    mid_table_list = []

    for line in lines:

        # parse input lines
        items = line.split(',')

        # each MID specification must have exactly 4 elements
        if len(items) != 4:
            print "\n Input file: ", in_file
            print " Line: ", line
            error("A MID specification must have exactly 4 elements")

        compound_name = items[0]
        rt = time_str_secs(items[1]) # convert to seconds
        diagnostic_ion = int(items[2])
        mdv_size = int(items[3])

        # set compound name, retention time, diagnostic ions and MDV size
        mid_table = MID_table(compound_name, rt, diagnostic_ion, mdv_size)

        # store MID table in MID table list
        mid_table_list.append(mid_table)

    return mid_table_list

def parse_data_defs(in_file):

    """
    @summary: Read data file names and return as a list

    @param in_file: The name of the file containing data file names
    @type in_file: StringType

    @return: The list of data file names
    @rtype: ListType

    @author: Milica Ng
    @author: Vladimir Likic
    """

    lines = file_lines(in_file, filter=True)

    data_files = []  
    for line in lines:
        data_files.append(line)

    return data_files

def write_mid_tables(mid_table_list, out_file):

    """
    @summary: Write MID tables, including any warnings, to out_file

    @param mid_table_list: List of MID tables 
    @type mid_table_list: ListType
    @param out_file: The name of the file used for writing
    @type out_file: StringType

    @return: None
    @rtype: NoneType

    @author: Milica Ng
    """

    for mid_table in mid_table_list:     
        mid_table.write(out_file)

def parse_mid_tables_file(in_file):

    """
    @summary: Read MID table file and return as a list of MID_table objects

    @param in_file: The name of the file containing ion definitions
    @type in_file: StringType

    @return: MID_table object list with compound names, 
        retention times, diagnostic ions and mdv size and values
    @rtype: ListType

    @author: Milica Ng
    """

    lines = file_lines(in_file, filter=True)
    mid_table_list = []
    mid_table = []

    for line in lines:

        # parse input lines
        items = line.split(',')

        # check if line is empty
        if len(items) != 1:

            # if line is not empty, check if it contains compound name
            if items[0] == 'compound':
                compound_name = items[1]
            # if line is not empty, check if it contains retention time
            if items[0] == 'rt =':
                rt = float(items[1])
            # if line is not empty, check if it contains diagnostic ion value
            if items[0] == 'ion':
                diagnostic_ion = int(items[1])
            # if line is not empty, check if it contains MID vector size
            if items[0] == 'mdv size =':
                mdv_size = int(items[1])
                # set compound name, retention time, diagnostic ion and MDV size
                mid_table = MID_table(compound_name, rt, diagnostic_ion, mdv_size)
            # if line is not empty, check if it contains file name 
            if items[0] == 'file name':
                file_name = items[1]
                mdv = []
                for i in range(2, mdv_size+2):
                    mdv.append(float(items[i]))
                # store MID values
                mid_table.set_values(mdv, file_name)
        else: 
            if not(mid_table == []):   
                # store MID table in MID table list
                mid_table_list.append(mid_table)
                # empty MID table
                mid_table = []

    #test_file = 'test.csv'
    #write_mid_tables(mid_table_list,test_file)

    return mid_table_list

def parse_ion_file(in_file, mid_table_list):

    """
    @summary: Read ion composition and return as part of MID_table objects list

    @param in_file: The name of the file containing ion compositions
    @type in_file: StringType
    @return: MID_table object list with compound names, 
        retention times, diagnostic ions and mdv size and values set
    @rtype: ListType

    @return: MID_table object list with ion compositions set 
    @rtype: ListType

    @author: Milica Ng
    """

    lines = file_lines(in_file, filter=True)

    for line in lines:

        # parse input lines
        items = line.split(',')

        # each MID specification must have exactly 4 elements
        if len(items) != 3:
            print "\n Input file: ", in_file
            print " Line: ", line
            error("Ion composition specification must have exactly 3 elements")

        compound_name = items[0]
        diagnostic_ion = int(items[1])
        ion_atom_composition = items[2]
        
        # split ion atom composition by spaces
        string = ion_atom_composition.split()

        for i in range(0,len(string)):

            word = string[i]
            element = ''
            num = ''

            for j in range(0,len(word)):
                char = word[j]
                if char >= 'a' and  char <= 'z':
                    element = element + char
                elif char >= 'A' and  char <= 'Z':
                    element = element + char.swapcase()
                elif char >= '0' and  char <= '9':
                    num = num + char 
                else:
                    error("Ion composition specification contains invalid character")  

            for mid_table in mid_table_list:
                if (mid_table.get_compound_name() == compound_name) and (mid_table.get_ion() == diagnostic_ion):
                    # print '\n element: ', element, '  num: ', num
                    mid_table.set_atoms(element, int(num))     

    return mid_table_list

