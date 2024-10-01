import numpy as np
import math
from stingray import Lightcurve, AveragedCrossspectrum

#calculating the rms and phase over modulation angle via the cross-spectrum of datasets

def real_im_G_calc(mod_min, mod_max,data_1, lc_ref,bin_length, seg_length, fmin, fmax,spur_sub=True):
    """ Process a single modulation angle range and return the averaged real and imaginary power. """
    # Filter data for the current modulation angle range
    data_phi=data['PHI']
    data_bin = data[(mod_min <= data_phi) & (data_phi <= mod_max)]
    time_bin = data_bin['TIME']
    
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

        Q_spur_1_norm=q_sum/I #define normalised spur stokes 
        U_spur_1_norm=u_sum/I
            
            
        spur_sub=((lc_countrate/mod_bin_number)* ((Q_spur_1_norm*np.cos((2*mod_midpoint)))+(U_spur_1_norm*np.sin(((2*mod_midpoint))))))*(bin_length) #assumes qsm,usm const over time
        spur_sub_counts=[spur_sub]*len(lc.time)
        lc_spur=Lightcurve(lc.time,spur_sub_counts)                                
        lc_counts_subtracted=lc.counts-spur_sub_counts #subtracting spur lc
        lc=Lightcurve(lc.time,lc_counts_subtracted)

    else:
        lc=lc


    # Create averaged cross spectrum
    cs = AveragedCrossspectrum.from_lightcurve(lc, lc_ref, seg_length, norm='frac')
    
    # Extract the real and imaginary parts of the power spectrum
    real_power = cs.power.real[(fmin <= cs.freq) & (cs.freq <= fmax)].mean()
    imag_power = cs.power.imag[(fmin <= cs.freq) & (cs.freq <= fmax)].mean()

    return real_power, imag_power



def rms_phase_calc(data_1, lc_ref,bin_length,seg_length, fmin, fmax,norm_factor,mod_bin_number,spur_sub=True):
    
    mod_min = np.radians(-90)
    mod_max = np.radians(90)
    aspace=np.linspace(mod_min,mod_max,mod_bin_number+1)
    mod_angle_list=[(aspace[i-1],aspace[i]) for i in range(1,len(aspace))]  #making a list of mod angle bins to select over
    mod_min_array=[i[0] for i in mod_angle_list]
    mod_max_array=[i[1] for i in mod_angle_list]
    av_mod=[np.mean(i) for i in mod_angle_list]
    av_mod_err=[(i[1]-i[0])/2 for i  in mod_angle_list]
    
    real_im_G_calc_VEC = np.vectorize(real_im_G_calc)
    #for i in mod_angle_list:
    #    mod_min_i, mod_max_i = i
        
        # Get the real and imaginary parts of the power spectrum for the current modulation range
    av_power_real, av_power_im = real_im_G_calc_VEC(mod_min_array, mod_max_array,data_1, lc_ref,bin_length, seg_length, fmin, fmax,spur_sub)
 
    # Convert real and imaginary power into complex G
    
    complex_G = [complex(a, b) for a, b in zip(av_power_array_real, av_power_array_im)]
    complex_G_arr.append(complex_G)
    # Calculate fractional RMS and phase
    frac_rms = [norm_factor * abs(i) for i in complex_G]
    phase = [math.atan2(i.imag, i.real) for i in complex_G]
    return frac_rms, phase


#Calculating Ingram 2019 Errorbars







