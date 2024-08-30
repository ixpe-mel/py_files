
import numpy as np
import time
import scipy.fft as fft
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.optimize import curve_fit
from astropy.io import ascii
from astropy.table import Table
import math
import stingray
import lightcurve
from stingray import Lightcurve
from stingray import Powerspectrum
from stingray import AveragedPowerspectrum
from stingray import Crossspectrum
from stingray.exceptions import StingrayError
from stingray import AveragedCrossspectrum
import scipy
from scipy import stats


#File 1 should be the combination of two DUs, usually 1 and 2

def Polarisation_variability_finder(file1,file2,output_file,gti,bin_length,seg_length,fmin,fmax,mod_bin_number):
    
    Pmin=51
    Pmax=200
    GTI=list(np.loadtxt(str(gti)))
  
    with fits.open(str(file1)) as hdu:
            data_1=hdu[1].data  #reading in first file of two detectors
    
    with fits.open(str(file2)) as hdu2:
            header=hdu2[1].header #reading in header 
            data_2=hdu2[1].data #reading in DU3
            TSTART=header['TSTART']
            TSTOP=header['TSTOP']

            
            # Data filtering
          
    #PI channel/energy index
    data_1=data_1[(Pmin<=data_1['PI']) & (data_1['PI']<= Pmax)]
    data_2=data_2[(Pmin<=data_2['PI']) & (data_2['PI']<= Pmax)]
    
    # Quality Factor cut
    
    data_1=data_1[data_1['QUAL']==1]
    data_2=data_2[data_2['QUAL']==1]

    TIME_1=data_1['TIME']
    TIME_2=data_2['TIME']
    
    
            # Normalisation Calculation

    #Lightcurves from file 1 

    lightcurve_1=Lightcurve.make_lightcurve(TIME_1,dt=bin_length,tstart=TSTART,gti=GTI)
    lightcurve_1.apply_gtis()
    
    
    # Lightcurve from file 2 (ditto)
    
    lightcurve_2=Lightcurve.make_lightcurve(TIME_2,dt=bin_length,tstart=TSTART,gti=GTI)
    lightcurve_2.apply_gtis()
    
    # Power spectrum from file 2
    
    ps_2=Powerspectrum.from_lightcurve(lightcurve_2,seg_length,norm='frac')
    ps_2_mean=ps_2.power[(fmin<=ps_2.freq) & (ps_2.freq<=fmax)].mean()

    #Cross spec for norm ie over all mod angle bins

    cs = AveragedCrossspectrum.from_lightcurve(lightcurve_1,lightcurve_2,seg_length,norm='frac')
    del lightcurve_1
    
    cs_real_mean=cs.power.real[(fmin<=cs.freq) & (cs.freq<=fmax)].mean()
    
    norm_factor=(np.sqrt((fmax-fmin))/np.sqrt(cs_real_mean))
    del cs_real_mean
    del cs
  
        
   
    # Making a list of mod angle bins to select over
    mod_minimum=np.radians(-90)
    mod_maximum=np.radians(90)
    
    aspace=np.linspace(mod_minimum,mod_maximum,mod_bin_number+1)
    mod_angle_list=[(aspace[i-1],aspace[i]) for i in range(1,len(aspace))]  #making a list of mod angle bins to select over
    
    
    
    # Calculating G for each modulation angle bin
    
    complex_G_array=[]
    ps_1_modbin_cut_mean_array=[]
    
    for i in mod_angle_list:
        mod_min=i[0]
        mod_max=i[1]
        mod_midpoint=np.mean(i)
        
       
        # Selecting photons within modulation angle bin range
        data_cut=data_1[(mod_min<=data_1['PHI']) & (data_1['PHI']<=mod_max)]
        
        # Creating lightcurve with selected photons
        lc=Lightcurve.make_lightcurve(data_cut['TIME'],dt=bin_length,gti=GTI)
        lc.apply_gtis()
        lc_countrate=lc.meanrate
        
        # Subtracting spurious polarisation
        
        q_spur_1=data_cut['QSP'] #defining spurious polarisation stokes parameters
        u_spur_1=data_cut['USP']

        q_sum=np.sum(q_spur_1)

        q_len=len(q_spur_1)

        Q_spur_1=np.sum(q_spur_1)/len(q_spur_1)
        U_spur_1=np.sum(u_spur_1)/len(u_spur_1)
        
        
        spur_sub=((lc_countrate/mod_bin_number)* ((Q_spur_1*np.cos((2*mod_midpoint)))+(U_spur_1*np.sin(((2*mod_midpoint))))))*(bin_length) #making spur lc
        spur_sub_counts=[spur_sub]*len(lc.time)
        lc_spur=Lightcurve(lc.time,spur_sub_counts)                                #REVISIT THEORY OF THIS
        lc_counts_subtracted=lc.counts-spur_sub_counts #subtracting spur lc
        lc=Lightcurve(lc.time,lc_counts_subtracted)
        
        # Creating Powerspectrum
        
        ps_1_modbin_cut=Powerspectrum.from_lightcurve(lc,seg_length,norm='frac')
        ps_1_modbin_cut_mean=ps_1_modbin_cut.power[(fmin<=ps_1_modbin_cut.freq) & (ps_1_modbin_cut.freq<=fmax)].mean()
        ps_1_modbin_cut_mean_array.append(ps_1_modbin_cut_mean)
        del ps_1_modbin_cut
        
        # Creating Cross spectrum
        
        cs=AveragedCrossspectrum.from_lightcurve(lc,lightcurve_2,seg_length,norm='frac')
        del lc
        cs_real_mean=cs.power.real[(fmin<=cs.freq) & (cs.freq<=fmax)].mean()
        cs_im_mean=cs.power.imag[(fmin<=cs.freq) & (cs.freq<=fmax)].mean()
        print(cs_real_mean)
        
        complex_G=tuple(zip(cs_real_mean, cs_im_mean))
        complex_G_array.append(complex_G)
        
            # Calculating results
    modulus_G=[np.abs(g) for g in complex_G_array]

    
            # Calculating Ingram 2019 errorbars
    n=len(cs.power.real[(fmin<=cs.freq) & (cs.freq<=fmax)])
    m=cs.m
    
    d_modulus_G=np.sqrt((1/(2*(n*m)))*((ps_1_modcut_mean_array*ps_2_mean)-((ps_2_mean/(cs_real_mean))*(modulus_G**2))))



    frac_rms=[i*norm_factor for i in modulus_G]
    phase=[(math.atan2(i))/(2*np.pi) for i in complex_G_array]
          
    frac_rms_err=norm_factor*d_modulus_G
    phase_err=d_modulus_G/modulus_G
    
    
    
    av_mod=[np.mean(i) for i in mod_angle_list]
    av_mod_err=[(i[1]-i[0])/2 for i  in mod_angle_list]
    
    results=np.array(tuple(zip(av_mod,av_mod_err,frac_rms,frac_rms_err,phase_phase_err)))
    np.savetxt(output_file,results)
          
    
        
          
    