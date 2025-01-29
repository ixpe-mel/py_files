import numpy as np
import math
from stingray import Lightcurve, AveragedCrossspectrum, Powerspectrum
from functools import partial
import warnings
import importlib
import matplotlib.pyplot as plt
warnings.filterwarnings('ignore')


def G_calc(mod_min, mod_max,data_1,lc_ref,GTI,bin_length, seg_length, fmin, fmax,mod_bin_number,spur_sub,norm):
    """ Process a single modulation angle range and return the averaged real and imaginary power. """
    # Filter data for the current modulation angle range

    av_mod=np.mean([mod_min,mod_max])
    data_phi=data_1['PHI']
    #print(len(data_phi))   
    data_bin = data_1[(mod_min <= data_phi) & (data_phi <= mod_max)]
    time_bin = data_bin['TIME']

    lc_1_ref=Lightcurve.make_lightcurve(data_1['TIME'], dt=bin_length, gti=GTI)
    lc_1_ref.apply_gtis()



       # Create the Subject band lightcurve
    lc = Lightcurve.make_lightcurve(time_bin, dt=bin_length, gti=GTI)
    lc.apply_gtis()
    lc_countrate=lc.meanrate

    lc_before=lc
   
    
    N=np.sum(lc_1_ref.counts)
    I=len(lc.counts)
    J=mod_bin_number

    mod_bin_width=(mod_max-mod_min)
    T=lc.time[-1] - lc.time[0]

    q_spur_1_TOT=data_1['QSP'] # all events
    u_spur_1_TOT=data_1['USP']

    Q_spur_1_TOT=np.sum(q_spur_1_TOT) # equivelent to big q
    U_spur_1_TOT=np.sum(u_spur_1_TOT) #equivalent to big u
    I_TOT=len(q_spur_1_TOT) # equivelent to I or N

    Q_spur_1_norm_TOT=np.sum(q_spur_1_TOT)/I_TOT #define normalised spur stokes 
    U_spur_1_norm_TOT=np.sum(u_spur_1_TOT)/I_TOT


    if spur_sub==True:
        #print('Subtracting spurious polarisation...')
        q_spur_1=data_bin['QSP'] # per event spurious stokes parameters for mod bin of interest
        u_spur_1=data_bin['USP']

        Q_spur_1=np.sum(q_spur_1) # equivelent to big q
        U_spur_1=np.sum(u_spur_1) #equivalent to big u
        I=len(q_spur_1) # equivelent to I or N

        Q_spur_1_norm=np.sum(q_spur_1)/I #define normalised spur stokes 
        U_spur_1_norm=np.sum(u_spur_1)/I
            

   #2nd VARIATION
        spur_sub=(bin_length/T)*(mod_bin_width/(np.pi))*((Q_spur_1_norm_TOT*np.cos(2*av_mod)) +(U_spur_1_norm_TOT*np.sin( 2*av_mod  )))
        #spur_sub=1e-5
   #3rd
        #spur_sub=(bin_length/T)*(mod_bin_width/np.pi)*((Q_spur_1_norm*np.cos(2*av_mod)) +(U_spur_1_norm*np.sin( 2*av_mod  )))


            
    #ORIGINAL 
    
        #spur_sub=((lc_countrate/mod_bin_number)* ((Q_spur_1_norm*np.cos((2*av_mod)))+(U_spur_1_norm*np.sin(((2*av_mod))))))*(bin_length) #assumes qsm,usm const over time
        spur_sub_counts=[spur_sub]*len(lc.time)
        #spur_sub_counts = [0.0]*len(lc.time)
        lc_spur=Lightcurve(lc.time,spur_sub_counts)     
        print('spur_sub counts',lc_spur.counts)                           
        lc_counts_subtracted=lc.counts-spur_sub_counts #subtracting spur lc
        lc_subject=Lightcurve(lc.time,lc_counts_subtracted)
        #plt.figure()
        #plt.title('lc spur')
        #plt.plot(lc_spur.time,lc_spur.counts,label='spur')
        #plt.legend()
        #plt.show()
        #do over time ? READ
        #are_same=lc_subject.counts==lc_before.counts
        #print('are_same',are_same)

    else:
        print('Not subtracting spurious polarisation...')
        lc_subject=lc
        lc_spur=0
    #print('lc before counts',lc_before.counts)
    #print('lc_subject counts',lc_subject.counts)
    #are_same=lc_subject.counts==lc_before.counts
    #print('are_same',are_same)

    #plt.figure()
    #plt.title('lc compare')
    #plt.plot(lc_before.time,lc_before.counts,label='before spur sub')
    #plt.plot(lc_subject.time,lc_subject.counts,label='after spur sub')
    #plt.show()

    # Create averaged cross spectrum
    cs = AveragedCrossspectrum.from_lightcurve(lc_subject, lc_ref, seg_length, norm='abs')
    #scale_factor=(2*(1/64))/(lc_subject.meanrate*lc_ref.meanrate*len(lc_subject.counts))

    cs_real_nu=cs.power.real[(fmin<=cs.freq) & (cs.freq<=fmax)]
    cs_imag_nu=cs.power.imag[(fmin<=cs.freq) & (cs.freq<=fmax)]
    cs_freq_nu=cs.freq[(fmin<=cs.freq) & (cs.freq<=fmax)]
    cs_sr_err_real=cs.power_err.real[(fmin<=cs.freq) & (cs.freq<=fmax)]
    cs_sr_err_imag=cs.power_err.imag[(fmin<=cs.freq) & (cs.freq<=fmax)]

    #PLotting real and im parts of power spectrum
    #print('mod_min',mod_min)
    #print('mod_max',mod_max)
    #plt.figure()
    #plt.title('real')
    #plt.xscale('log')
    #plt.yscale('log')
    #plt.errorbar(cs_freq_nu,cs_real_nu,yerr=cs_sr_err_real,fmt='o',label='Real')
    #plt.show()
    #plt.figure()
    #plt.title('imag')
    #plt.xscale('log')
    #plt.yscale('log')
    #plt.errorbar(cs_freq_nu,cs_imag_nu,yerr=cs_sr_err_imag,fmt='o',label='Imaginary')
    #plt.show()

    cs_sr_err_real=np.sqrt(np.sum(cs_sr_err_real**2))/len(cs_sr_err_real)
    #print('cs_sr_err_real',cs_sr_err_real)
    cs_sr_err_imag=np.sqrt(np.sum(cs_sr_err_imag**2))/len(cs_sr_err_imag)
    #print('cs_sr_err_imag',cs_sr_err_imag)
    # Extract the real and imaginary parts of the power spectrum
    G_real = cs.power.real[(fmin <= cs.freq) & (cs.freq <= fmax)].mean()
    #print('G_real',G_real)
    G_im =cs.power.imag[(fmin <= cs.freq) & (cs.freq <= fmax)].mean()

    #f_mid=(fmin+fmax)/2

    #plt.figure()
    #plt.title('Greal')
    #plt.xscale('log')
    #plt.yscale('log')
    #plt.errorbar(f_mid,G_real,yerr=cs_sr_err_real,fmt='o',label='Real')
    #plt.show()
    #plt.figure()
    #plt.title('Gimag')
    #plt.xscale('log')
    #plt.yscale('log')
    #plt.errorbar(f_mid,G_im,yerr=cs_sr_err_imag,fmt='o',label='Imaginary')
    #plt.show()


    n=len(cs.power.real[(fmin<=cs.freq) & (cs.freq<=fmax)]) #number of datapoints of power spectrum in frequency range
    m=cs.m #the number of cross spectra averaged together
    
    return G_real, G_im, n, m, lc_subject,lc_spur,cs