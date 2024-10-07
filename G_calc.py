import numpy as np
import math
from stingray import Lightcurve, AveragedCrossspectrum, Powerspectrum
from functools import partial
import warnings
warnings.filterwarnings('ignore')


def G_calc(mod_min, mod_max,data_1,data_2,lc_ref,GTI,bin_length, seg_length, fmin, fmax,mod_bin_number,spur_sub=True,coherence_corrector=True):
    """ Process a single modulation angle range and return the averaged real and imaginary power. """
    # Filter data for the current modulation angle range
    #print(mod_min,mod_max)
    av_mod=np.mean([mod_min,mod_max])
    data_phi=data_1['PHI']
    #print(len(data_phi))   
    data_bin = data_1[(mod_min <= data_phi) & (data_phi <= mod_max)]
    time_bin = data_bin['TIME']
    #print(time_bin)
    #print(GTI)
       # Create the lightcurve
    lc = Lightcurve.make_lightcurve(time_bin, dt=bin_length, gti=GTI)
    lc.apply_gtis()
    lc_countrate=lc.meanrate

   

    if spur_sub==True:
        print('Subtracting spurious polarisation...')
        q_spur_1=data_bin['QSP'] # per event spurious stokes parameters for mod bin of interest
        u_spur_1=data_bin['USP']

        Q_spur_1=np.sum(q_spur_1) # equivelent to big q
        U_spur_1=np.sum(u_spur_1) #equivalent to big u
        I=len(q_spur_1) # equivelent to I or N

        Q_spur_1_norm=np.sum(q_spur_1)/I #define normalised spur stokes 
        U_spur_1_norm=np.sum(u_spur_1)/I
            
            
        spur_sub=((lc_countrate/mod_bin_number)* ((Q_spur_1_norm*np.cos((2*av_mod)))+(U_spur_1_norm*np.sin(((2*av_mod))))))*(bin_length) #assumes qsm,usm const over time
        spur_sub_counts=[spur_sub]*len(lc.time)
        lc_spur=Lightcurve(lc.time,spur_sub_counts)                                
        lc_counts_subtracted=lc.counts-spur_sub_counts #subtracting spur lc
        lc_subject=Lightcurve(lc.time,lc_counts_subtracted)

    else:
        print('Not subtracting spurious polarisation...')
        lc_subject=lc


    # Create averaged cross spectrum
    cs = AveragedCrossspectrum.from_lightcurve(lc_subject, lc_ref, seg_length, norm='frac')
    
    # Extract the real and imaginary parts of the power spectrum
    G_real = cs.power.real[(fmin <= cs.freq) & (cs.freq <= fmax)].mean()
    G_im = cs.power.imag[(fmin <= cs.freq) & (cs.freq <= fmax)].mean()

    n=len(cs.power.real[(fmin<=cs.freq) & (cs.freq<=fmax)]) #number of datapoints of power spectrum in frequency range
    m=cs.m #the number of cross spectra averaged together
    return G_real, G_im, n, m,lc_subject
    