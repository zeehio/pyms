
import sys
sys.path.append('/x/PyMS/')

import string

from pyms.MSlib.Class import msrecord
from pyms.MSlib.Class import mslib
from pyms.Utils.IO import file_lines
from pyms.Utils.IO import dump_object
from pyms.Utils.IO import load_object


def load_nist(file_name):    

    __CMPD_NAME_KEYWORD = "##CAS NAME"
    __CAS_REG_KEYWORD = "##CAS REGISTRY NO"
    __NUM_PEAKS_KEYWORD = "##NPOINTS"
    __RECORD_START_KEYWORD = "##TITLE"

    nist_lines = file_lines(file_name)
    
    collect_ms_flag = False
    cmpd_counter = 0
    ms_lib = mslib()

    for line in nist_lines:        
        fields = string.split(line, "=")

        if len(fields)>0 and fields[0].upper() == __RECORD_START_KEYWORD:
            collect_ms_flag = False
            
        elif len(fields)>0 and fields[0].upper() == __CMPD_NAME_KEYWORD:
            ms_record = msrecord()
            cmpd_counter = cmpd_counter+1            

            if(cmpd_counter>1):
                ms_record.set(cmpd_name, cmpd_regno, mi_dict)
                ms_lib.addrecord(ms_record)

            keyword_value = fields[1]
            cmpd_name = keyword_value.strip()

        elif len(fields)>0 and fields[0].upper() == __CAS_REG_KEYWORD:
            keyword_value = fields[1]
            cmpd_regno = keyword_value.strip()

        elif len(fields)>0 and fields[0].upper() == __NUM_PEAKS_KEYWORD:
            keyword_value = fields[1]
            num_peaks_str = keyword_value.strip()
            num_peaks = int(num_peaks_str)
            collect_ms_flag = True
            mi_dict = {}

        elif collect_ms_flag:
            if(fields[0][0]!='#'):
                mi = string.split(line)
                mass = int(mi[0])
                intensity = int(mi[1])
                mi_dict[mass] = intensity

    ms_record.set(cmpd_name, cmpd_regno, mi_dict)
    ms_lib.addrecord(ms_record)
       
    return ms_lib


def write_ms_lib(ms_lib, filename):
    dump_object(ms_lib, filename)
                
                
def load_ms_lib(filename):
    print "Loading ms library object..."
    lib = load_object(filename)
    return lib

def read_ms_lib(ms_lib):
    counter = 0
    for x in ms_lib.msrecord_list:
        counter = counter + 1
        print counter, x.name    
    
    