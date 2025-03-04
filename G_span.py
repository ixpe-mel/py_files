#Vectorizing G_calc to produce evenly spaced G from a list of modulation angles
import G_calc as Gc
#import dG_calc as dG_calc
import numpy as np
from functools import partial
import importlib
importlib.reload(Gc)
def G_span(mod_bin_number,data_1,lc_ref,GTI,bin_length, seg_length, fmin, fmax,spur_sub,norm='frac'):

    #create modulation angle bins
    mod_min_global = np.radians(-90)
    mod_max_global = np.radians(90)
    aspace = np.linspace(mod_min_global, mod_max_global, mod_bin_number + 1)

    mod_min_array = aspace[:-1]
    mod_max_array = aspace[1:]
    av_mod = (mod_min_array + mod_max_array) / 2
    av_mod_err = (mod_max_array - mod_min_array) / 2
    #print('mod_min_array',mod_min_array)
    #print('mod_max_array',mod_max_array)
    #Make partial function
    partial_G_calc = partial(Gc.G_calc, data_1=data_1,lc_ref=lc_ref,GTI=GTI,bin_length=bin_length, seg_length=seg_length, fmin=fmin, fmax=fmax,mod_bin_number=mod_bin_number,spur_sub=spur_sub,norm=norm)
    
    #vectorize real/im calculator
    partial_G_calc_VEC = np.vectorize(partial_G_calc,otypes=[object])

    # Get the real and imaginary parts of the power spectrum for the current modulation range
    result= partial_G_calc_VEC(mod_min_array, mod_max_array)#,data_1, lc_ref,GTI,bin_length, seg_length, fmin, fmax,spur_sub)
    G_real, G_im, n, m, lc_subject,lc_spur,cs,spur_sub_norm= zip(*result)
    #cs_err_sr_real,cs_err_sr_im ,cs= zip(*result)
    return G_real, G_im, n, m,lc_subject,lc_spur,cs,spur_sub_norm
#,lc_sub
#,cs_err_sr_real,cs_err_sr_im,cs