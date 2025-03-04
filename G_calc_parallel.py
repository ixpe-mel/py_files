#Does the same thing as G_calc.py but instead of reading in a simulated data file,
#it reads in the separated requirements so we dont need to make
#a simulated fits file.
import numpy as np
import math
from stingray import Lightcurve, AveragedCrossspectrum, Powerspectrum
from functools import partial
import warnings
import importlib

warnings.filterwarnings('ignore')
def G_calc(mod_min, mod_max,data_1_times,data_1_phi,data_1_qsp,data_1_usp,
           lc_ref,GTI,bin_length, seg_length, fmin, fmax,mod_bin_number,spur_sub,coherence_corrector):
    
    
    
    """ Process a single modulation angle range and return the averaged real and imaginary power. """
    # Filter data for the current modulation angle range
    #print('start')
    #print(np.degrees(mod_min),np.degrees(mod_max))
    data_1=[data_1_times,data_1_qsp,data_1_usp]
    lc_1_ref=Lightcurve.make_lightcurve(data_1_times, dt=bin_length, gti=GTI)
    lc_1_ref.apply_gtis()


    #print('times',data_1[0])
    #print('qsp',data_1[1])
    #print('usp',data_1[2])
    av_mod=np.mean([mod_min,mod_max])
    #data_phi=data_1['PHI']
    #print(len(data_phi))   
    data_bin = [arr[(mod_min <= data_1_phi) & (data_1_phi <= mod_max)] for arr in data_1]
    time_bin = data_bin[0]
    
   

    #print('time bin',time_bin)
    #print(GTI)
       # Create the lightcurve
    lc = Lightcurve.make_lightcurve(time_bin, dt=bin_length, gti=GTI)
    lc.apply_gtis()
    lc_countrate=lc.meanrate
    mod_bin_width=(mod_max-mod_min)
    T=lc.time[-1] - lc.time[0]

    N=np.sum(lc_1_ref.counts)
    I=len(lc.counts)
    J=mod_bin_number

    q_spur_1_TOT=data_1_qsp # all events
    u_spur_1_TOT=data_1_usp

    Q_spur_1_TOT=np.sum(q_spur_1_TOT) # equivelent to big q
    U_spur_1_TOT=np.sum(u_spur_1_TOT) #equivalent to big u
    I_TOT=len(q_spur_1_TOT) # equivelent to I or N

    Q_spur_1_norm_TOT=np.sum(q_spur_1_TOT)/I_TOT #define normalised spur stokes 
    U_spur_1_norm_TOT=np.sum(u_spur_1_TOT)/I_TOT

   

    if spur_sub==True:
        print('Subtracting spurious polarisation...')
        q_spur_1=data_bin[1] # per event spurious stokes parameters for mod bin of interest
        u_spur_1=data_bin[2]

        Q_spur_1=np.sum(q_spur_1) # equivelent to big q
        U_spur_1=np.sum(u_spur_1) #equivalent to big u
        I=len(q_spur_1) # equivelent to I or N

        Q_spur_1_norm=np.sum(q_spur_1)/I #define normalised spur stokes 
        U_spur_1_norm=np.sum(u_spur_1)/I
            
        #2nd VARIATION
        spur_sub=(bin_length/T)*(mod_bin_width/np.pi)*((Q_spur_1_norm_TOT*np.cos(2*av_mod)) +(U_spur_1_norm_TOT*np.sin( 2*av_mod  )))
        
        
        
        # ORIGINAL spur_sub=((lc_countrate/mod_bin_number)* ((Q_spur_1_norm*np.cos((2*av_mod)))+(U_spur_1_norm*np.sin(((2*av_mod))))))*(bin_length) #assumes qsm,usm const over time
        spur_sub_counts=[spur_sub]*len(lc.time)
        lc_spur=Lightcurve(lc.time,spur_sub_counts)                                
        lc_counts_subtracted=lc.counts-spur_sub_counts #subtracting spur lc
        lc_subject_spur_subtracted=Lightcurve(lc.time,lc_counts_subtracted)

        lc_subject=lc
        cs=AveragedCrossspectrum.from_lightcurve(lc,lc_ref,seg_length,norm='frac')
        normalisation=lc.meanrate/(lc.meanrate-lc_spur.meanrate)

   
        #do over time ? READ

    else:
        #print('Not subtracting spurious polarisation...')
        lc_subject=lc
        lc_spur=0
        normalisation=1
        cs=AveragedCrossspectrum.from_lightcurve(lc,lc_ref,seg_length,norm='frac')
    
    

    #print('lc_subject',lc_subject)
    #print('lc_ref',lc_ref)
    # Create averaged cross spectrum
    #cs = AveragedCrossspectrum.from_lightcurve(lc_subject, lc_ref, seg_length, norm='frac')
    #print('real',cs.power.real)
    #print('imag',cs.power.imag)
    #print('freq',cs.freq)
    #print('fmin',fmin)
    #print('fmax',fmax)
    index=[(fmin <= cs.freq) & (cs.freq <= fmax)]
    #print('index',index)


    #print("cs.freq:", cs.freq)
    #print("cs.power.real:", cs.power.real)
    #print("Selected freq range:", cs.freq[(fmin <= cs.freq) & (cs.freq <= fmax)])
    #print("Mean of selected power values:", cs.power.real[(fmin <= cs.freq) & (cs.freq <= fmax)])

    #np.savetxt('/home/c2032014/Inject_signal/freq_test.txt',cs.freq)
    G_real = cs.power.real[(fmin <= cs.freq) & (cs.freq <= fmax)].mean()
    G_real=normalisation*G_real
    # Extract the real and imaginary parts of the power spectrum
   
    #print('G_real',G_real)
    G_im = cs.power.imag[(fmin <= cs.freq) & (cs.freq <= fmax)].mean()
    G_im=normalisation*G_im

    n=len(cs.power.real[(fmin<=cs.freq) & (cs.freq<=fmax)]) #number of datapoints of power spectrum in frequency range
    m=cs.m #the number of cross spectra averaged together
    return G_real, G_im, n, m, lc_subject
    