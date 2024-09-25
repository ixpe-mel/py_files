
# Make a power spectrum from an IXPE fits file and GTI list. 

# Bin length= length of time bin used to create lightcurve
# Seg length = Length of segments averaged over in the power spectrum

import numpy as np
import scipy.fft as fft
from astropy.io import fits
import math
from stingray import AveragedPowerspectrum
from stingray import Lightcurve
import matplotlib.pyplot as plt
def powerspectrum(data_file,gti,bin_length,seg_length):

    Pmin=51
    Pmax=200
    GTI=list(np.loadtxt(str(gti))) #read in GTIs
    
    
    with fits.open(str(data_file)) as hdu:
        data=hdu[1].data
        
        #PI channel cut
        data_=data[(Pmin<=data['PI']) & (data['PI']<= Pmax)]
  
        TIME=data.field('TIME')
       
        # Create Lightcurve

        lightcurve_12=Lightcurve.make_lightcurve(TIME,dt=bin_length,gti=GTI)
        lightcurve_12.apply_gtis()

        
        #Create Power spectrum
        
        ps=AveragedPowerspectrum.from_lightcurve(lightcurve_12,seg_length,norm='frac')
        
        ps_err=ps.power_err
       

        fig, ax1 = plt.subplots(1,1,figsize=(9,6))
        ax1.plot(ps.freq, ps.power.real, color='green')
        #ax1.errorbar(ps.freq, ps.power.real, yerr=ps.power_err)        
            
        #ax1.plot(ps_log.freq, ps_log.power, color='red',label='log rebin')
        ax1.set_xlabel("Frequency (Hz)")
        ax1.set_title('Power spec')
        ax1.set_ylabel("Real Power")
        ax1.set_yscale('log')
        ax1.set_xscale('log')
        ax1.tick_params(axis='x', labelsize=16)
        #ax1.set_xlim(0,0.004938271604+1)
        #ax1.set_ylim(0.0027,0.15)
        ax1.tick_params(axis='y', labelsize=16)
        ax1.tick_params(which='major', width=1.5, length=7)
        ax1.tick_params(which='minor', width=1.5, length=4)
        for axis in ['top', 'bottom', 'left', 'right']:
            
    
            ax1.spines[axis].set_linewidth(1.5)
            plt.show()

            
        result=np.array(tuple(zip(ps.freq,ps.power,ps.power*ps.freq,ps_err)))
        return ps.freq, ps.power, ps.power*ps.freq, ps_err
