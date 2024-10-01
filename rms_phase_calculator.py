import numpy as np
import math
from stingray import Lightcurve, AveragedCrossspectrum

#calculating the rms and phase over modulation angle via the cross-spectrum of datasets

# Assuming these variables are defined elsewhere in your code
# mod_angle_list, data_12, mod_angle_array_12, bin_length, GTI, lc_ref, seg_length, fmin, fmax, norm_factor

def process_mod_angle(mod_min, mod_max, mod_angle_array, data, lc_ref, seg_length, fmin, fmax):
    """ Process a single modulation angle range and return the averaged real and imaginary power. """
    # Filter data for the current modulation angle range
    data_bin = data[(mod_min <= mod_angle_array) & (mod_angle_array <= mod_max)]
    time_bin = data_bin['TIME']
    
    # Create the lightcurve
    lc = Lightcurve.make_lightcurve(time_bin, dt=bin_length, gti=GTI)
    lc.apply_gtis()
    
    # Create averaged cross spectrum
    cs = AveragedCrossspectrum.from_lightcurve(lc, lc_ref, seg_length, norm='frac')
    
    # Extract the real and imaginary parts of the power spectrum
    real_power = cs.power.real[(fmin <= cs.freq) & (cs.freq <= fmax)].mean()
    imag_power = cs.power.imag[(fmin <= cs.freq) & (cs.freq <= fmax)].mean()

    return real_power, imag_power



def rms_phase_calc(mod_min, mod_max, mod_angle_array_12, data_12, lc_ref, seg_length, fmin, fmax,norm_factor):
    for i in mod_angle_list:
        mod_min, mod_max = i
        
        # Calculate the average modulation value and append to the list
        mod = np.mean(i)
        mod_arr.append(mod)

        # Get the real and imaginary parts of the power spectrum for the current modulation range
        av_power_real, av_power_im = process_mod_angle(mod_min, mod_max, mod_angle_array_12, data_12, lc_ref, seg_length, fmin, fmax)
        
        # Append the results to the respective arrays
        av_power_array_real.append(av_power_real)
        av_power_array_im.append(av_power_im)

    # Convert real and imaginary power into complex numbers
    
    complex_G = [complex(a, b) for a, b in zip(av_power_array_real, av_power_array_im)]

    # Calculate fractional RMS and phase
    frac_rms = [norm_factor * abs(i) for i in complex_G]
    phase = [math.atan2(i.imag, i.real) for i in complex_G]

    return mod_arr, frac_rms, phase