"""
Provides helper functions for MIDs processing
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
from pyms.MID.Class import MIDS
from pyms.Utils.IO import file_lines

def parse_input(in_file):

    lines = file_lines(in_file)
    mids_list = []

    for line in lines:

        # parse input file
        items =line.split(',')
        name = str(items[0])
        rt = float(items[1])*60 # convert to seconds
        ion = int(items[2])
        mid_size = int(items[3])

        # set compound name, retention time, diagnostic ions and MID size
        mids = MIDS(name, rt, ion, mid_size)

        # store mids in mids_list
        mids_list.append(mids)

    return mids_list

