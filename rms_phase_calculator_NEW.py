import numpy as np
import math
from stingray import Lightcurve, AveragedCrossspectrum, Powerspectrum
from functools import partial
import Ingram_2019_errors as I_19errs
import warnings
warnings.filterwarnings('ignore')
#calculating the rms and phase over modulation angle via the cross-spectrum of datasets

def G_calc(mod_min, mod_max,av_mod,data_1,data_2,lc_ref,ps_2_mean,GTI,bin_length, seg_length, fmin, fmax,mod_bin_number,spur_sub=True,coherence_corrector=True):
    """ Process a single modulation angle range and return the averaged real and imaginary power. """
    # Filter data for the current modulation angle range
    #print(mod_min,mod_max)
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
        lc_subject=lc


    # Create averaged cross spectrum
    cs = AveragedCrossspectrum.from_lightcurve(lc_subject, lc_ref, seg_length, norm='frac')
    
    # Extract the real and imaginary parts of the power spectrum
    G_real = cs.power.real[(fmin <= cs.freq) & (cs.freq <= fmax)].mean()
    G_im = cs.power.imag[(fmin <= cs.freq) & (cs.freq <= fmax)].mean()

    n=len(cs.power.real[(fmin<=cs.freq) & (cs.freq<=fmax)]) #number of datapoints of power spectrum in frequency range
    m=cs.m #the number of cross spectra averaged together
    return G_real, G_im, n, m
    

def dG_calc(mod_min,mod_max,G_real,G_im,lc_subject,lc_ref,data_2,ps_2_mean,n,m,fmin,fmax,seg_length,bin_length,GTI,coherence_corrector)

    ps_2=Powerspectrum.from_lightcurve(lc_ref,seg_length,norm='frac')
    ps_2_mean=ps_2.power[(fmin<=ps_2.freq) & (ps_2.freq<=fmax)].mean()

    if coherence_corrector==True:
        print('Applying coherence correction...')
        dG=I_19errs.ingr_2019_errs_cc(mod_min,mod_max,G_real,G_im,lc_subject,data_2,ps_2_mean,n,m,fmin,fmax,seg_length,bin_length,GTI)
    else:
        print('Applying no coherence correction...')
        dG=I_19errs.ing_2019_errs(G_real,G_im,lc_subject,ps_2_mean,n,m,fmin,fmax,seg_length,bin_length,GTI)
    
    return dG



def rms_phase_anderr_calc(data_1,data_2,gti, lc_ref,bin_length,seg_length, fmin, fmax,norm_factor,mod_bin_number,spur_sub=True):
    
    #load in GTI
    GTI=list(np.loadtxt(str(gti)))

    #create modulation angle bins
    mod_min = np.radians(-90)
    mod_max = np.radians(90)
    aspace=np.linspace(mod_min,mod_max,mod_bin_number+1)
    mod_angle_list=[(aspace[i-1],aspace[i]) for i in range(1,len(aspace))]  #making a list of mod angle bins to select over
    mod_min_array=[i[0] for i in mod_angle_list]
    mod_max_array=[i[1] for i in mod_angle_list]
    av_mod=[np.mean(i) for i in mod_angle_list]
    av_mod_err=[(i[1]-i[0])/2 for i  in mod_angle_list]


    
    #Make partial function
    partial_G_calc = partial(G_dG_calc, data_1=data_1 ,data_2=data_2,ps_2_mean=ps_2_mean,lc_ref=lc_ref,GTI=GTI,bin_length=bin_length, seg_length=seg_length, fmin=fmin, fmax=fmax,mod_bin_number=mod_bin_number,spur_sub=spur_sub)
    partial_dG_calc=partial(dG_calc,data_2=data_2,ps_2_mean=ps_2_mean,n=n,m=m,fmin=fmin,fmax=fmax,seg_length=seg_length,bin_length=bin_length,GTI=GTI,coherence_corrector=coherence_corrector)
    #vectorize real/im calculator
    partial_G_calc_VEC = np.vectorize(partial_real_im_G_calc)
    partial_dG_calc_VEC=np.vectorize(partial_dG_calc)
        
    # Get the real and imaginary parts of the power spectrum for the current modulation range
    G_real, G_im = partial_G_calc_VEC(mod_min_array, mod_max_array,av_mod)#,data_1, lc_ref,GTI,bin_length, seg_length, fmin, fmax,spur_sub)
    dG=dG_calc_VEC(mod_min_array,mod_max_array,G_real,G_im,lc_subject,data_2,ps_2_mean,n,m,fmin,fmax,seg_length,bin_length,GTI,coherence_corrector)
    # Convert real and imaginary power into complex G
    
    complex_G = [complex(a, b) for a, b in zip(G_real, G_im)]
    # Calculate fractional RMS and phase
    frac_rms = [norm_factor * abs(i) for i in complex_G]
    phase = [math.atan2(i.imag, i.real) for i in complex_G]
    
    frac_rms_err=norm_factor*dG
    abs_vec=np.vectorize(abs)
    phase_err=dG/abs_vec(complex_G)

    return av_mod,av_mod_err,frac_rms,frac_rms_err ,phase,phase_err










