import matplotlib.pyplot as plt

import sys
sys.path.append('/x/PyMS/')

from pyms.Utils.Error import error
from pyms.GCMS.Class import IonChromatogram 


def plot_ics(ics, legend=None, label=None):
    """
    @summary: Plots an Ion Chromatogram or List of same
    
    @param ics: The ion chromatogram or list of same
    @type ics: pyms.GCMS.Class.IonChromatogram
    
    @param legend: Information for plot legend
    @type legend: String or list of strings
    
    @param label: A label for the plot
    @type label: String Type
    """
    
    # Container to store plots
    ic_plots = []
    		
    # color dictionary for plotting of ics; blue reserved
    # for TIC
    col_ic = {0:'b', 1:'r', 2:'g', 3:'k', 4:'y', 5:'m', 6:'c'}
   
    			
    #Plotting Variables
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    if not isinstance(ics, list):
        if isinstance(ics, IonChromatogram):
		ics = [ics]

        else:
            error("ics argument must be an IonChromatogram\
            or a list of Ion Chromatograms")
    # TODO: take care of case where one element of ics is
    # not an IonChromatogram
    
    if not isinstance(legend, list):
        legend = [legend]
		
	
    intensity_list = []
    time_list = ics[0].get_time_list()
			
	
    for i in range(len(ics)):
        intensity_list.append(ics[i].get_intensity_array())
    	
    for i in range(len(ics)):	
    	
        ic_plots.append(plt.plot(time_list, \
        intensity_list[i], col_ic[i]\
        , label = legend[i]))
        
        
        
    if label != None :
	t = ax.set_title(label)
	
    l = ax.legend()
			
    fig.canvas.draw
    plt.show()
    
