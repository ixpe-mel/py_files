# Finding the spin period of a pulsar by searching for the spin that gives the highest standard deviation between points of a lightcurve folded over spin period

import numpy as np
from astropy.io import fits
import math
import sys
sys.path.append('/home/c2032014/PhD/py_files')
import spin_phase_calculator as phase
import matplotlib.pyplot as plt

def find_pulsar_spin_period(fits_file,v_min,v_max,v_step):
    #create array of spin frequencies to look through
    
    with fits.open(str(fits_file)) as hdu:
         data=hdu[1].data
         header=hdu[1].header
         TSTART=header['TSTART']
         TIMES=data['TIME']  #Read in fits file
         times_sort=np.array(sorted(TIMES))


    floor_vec=np.vectorize(math.floor) #vectorise floor function

    
    t_0=TSTART #define arbitrary t_0 for phase calculation

    
    v_array=np.arange(v_min,v_max,v_step)
    
                   
                   
    phase_sd_array=[] 
                   
                   
    for v_spin in v_array:
    
    
        phase_i=phase.phase(v_spin,times_sort,t_0)
                   
        phase_i_sd=np.std(phase_i)
        phase_sd_array.append(phase_i_sd)
                   
    plt.figure()
    plt.plot(v_array,phase_sd_array,'.')
    plt.title('Standard deviation search for matching Neutron Star spin frequency')
    plt.xlabel('Spin frequency v')
    plt.ylabel('Standard deviation')
    
    max_index_sd=phase_sd_array.index(max(phase_sd_array))
    
    spin_freq=v_array[max_index_sd]
    
    print('Spin frequency:', spin_freq)
    
    
    
    
    
    
    
    
def find_pulsar_spin_period_higher_order(fits_file,v_min,v_max,v_step):
    #create array of spin frequencies to look through
    
    with fits.open(str(fits_file)) as hdu:
         data=hdu[1].data
         header=hdu[1].header
         TSTART=header['TSTART']
         TIMES=data['TIME']  #Read in fits file
         times_sort=np.array(sorted(TIMES))


    floor_vec=np.vectorize(math.floor) #vectorise floor function

    
    t_0=TSTART #define arbitrary t_0 for phase calculation

    
    v_array=np.arange(v_min,v_max,v_step)
    
                   
                   
    phase_sd_array=[] 
                   
                   
    for v_spin in v_array:
    
    
        phase_i=phase.phase(v_spin,times_sort,t_0)
                   
        phase_i_sd=np.std(phase_i)
        phase_sd_array.append(phase_i_sd)
                   
    plt.figure()
    plt.plot(v_array,phase_sd_array,'.')
    plt.title('Standard deviation search for matching Neutron Star spin frequency')
    plt.xlabel('Spin frequency v')
    plt.ylabel('Standard deviation')
    
    max_index_sd=phase_sd_array.index(max(phase_sd_array))
    
    spin_freq=v_array[max_index_sd]
    
    print('Spin frequency:', spin_freq)
