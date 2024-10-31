
#G span compatible with parallel processing
#Vectorizing G_calc to produce evenly spaced G from a list of modulation angles
import G_calc_parallel as Gcp
#import dG_calc as dG_calc
import numpy as np
from functools import partial
import importlib
importlib.reload(Gc)
def G_span(mod_min,mod_max,mod_bin_number,data_1_times,data_1_phi,data_1_qsp,data_1_usp,lc_ref,GTI,bin_length, seg_length, fmin, fmax,spur_sub=False,coherence_corrector=False):

    #create modulation angle bins
 
    aspace = np.linspace(mod_min, mod_max, mod_bin_number + 1)
    mod_min_array = aspace[:-1]
    mod_max_array = aspace[1:]
    av_mod = (mod_min_array + mod_max_array) / 2
    av_mod_err = (mod_max_array - mod_min_array) / 2

    #Make partial function
    partial_G_calc = partial(Gcp.G_calc, 
                             data_1_times=data_1_times,data_1_phi=data_1_phi,data_1_qsp=data_1_qsp,data_1_usp=data_1_usp,
                             lc_ref=lc_ref,GTI=GTI,bin_length=bin_length, seg_length=seg_length, fmin=fmin, fmax=fmax,mod_bin_number=mod_bin_number,spur_sub=spur_sub,coherence_corrector=coherence_corrector)
    
    #vectorize real/im calculator
    partial_G_calc_VEC = np.vectorize(partial_G_calc,otypes=[object])

    # Get the real and imaginary parts of the power spectrum for the current modulation range
    result= partial_G_calc_VEC(mod_min_array, mod_max_array)#,data_1, lc_ref,GTI,bin_length, seg_length, fmin, fmax,spur_sub)
    G_real, G_im, n, m, lc_subject = zip(*result)
    return G_real, G_im, n, m,lc_subject