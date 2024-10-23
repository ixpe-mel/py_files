#Calculating dG over a range of modulation angles, similarlly to G_span
import numpy as np
import dG_calc as dGc
from functools import partial
import importlib
importlib.reload(dGc)
def dG_span(mod_bin_number,G_real,G_im,lc_subject,lc_1_ref,lc_2_ref,data_2,n,m,fmin,fmax,seg_length,bin_length,GTI,coherence_corrector=True,spurious_sub=True):

    #create modulation angle bins
    mod_min_global = np.radians(-90)
    mod_max_global = np.radians(90)
    aspace = np.linspace(mod_min_global, mod_max_global, mod_bin_number + 1)

    mod_min_array = aspace[:-1]
    mod_max_array = aspace[1:]
    av_mod = (mod_min_array + mod_max_array) / 2
    av_mod_err = (mod_max_array - mod_min_array) / 2

    if coherence_corrector==True and spurious_sub==True:
        print('Applying coherence correction and spurious sub...')
        partial_dG_calc=partial(dGc.dG_calc, G_real=G_real,
               G_im=G_im,lc_1_ref=lc_1_ref,lc_2_ref=lc_2_ref,data_2=data_2,n=n,m=m,fmin=fmin,fmax=fmax,seg_length=seg_length,bin_length=bin_length,mod_bin_number=mod_bin_number,GTI=GTI)
        partial_dG_calc_VEC = np.vectorize(partial_dG_calc,otypes=[object])

        result= partial_dG_calc_VEC(mod_min_array, mod_max_array,lc_subject)
    
    
    elif coherence_corrector==True and spurious_sub==False:
        print('Applying coherence correction BUT NO SPUR SUB...')
        partial_dG_calc=partial(dGc.dG_calc, G_real=G_real,
               G_im=G_im,lc_1_ref=lc_1_ref,lc_2_ref=lc_2_ref,data_2=data_2,n=n,m=m,fmin=fmin,fmax=fmax,seg_length=seg_length,bin_length=bin_length,GTI=GTI,coherence_corrector=False)
        partial_dG_calc_VEC = np.vectorize(partial_dG_calc,otypes=[object])
        print(partial_dG_calc_VEC)
        result = partial_dG_calc_VEC(mod_min_array, mod_max_array,lc_subject)
        print(result)
    
    
    elif coherence_corrector==False and spurious_sub==True:

        partial_dG_calc=partial(dGc.dG_calc, G_real=G_real,
               G_im=G_im,lc_1_ref=lc_1_ref,lc_2_ref=lc_2_ref,data_2=data_2,n=n,m=m,fmin=fmin,fmax=fmax,seg_length=seg_length,bin_length=bin_length,mod_bin_number=mod_bin_number,GTI=GTI,spurious_sub=False)
        partial_dG_calc_VEC = np.vectorize(partial_dG_calc,otypes=[object])
        result= partial_dG_calc_VEC(mod_min_array, mod_max_array,lc_subject)
    
    
    #
    else:    
        partial_dG_calc=partial(dGc.dG_calc, G_real=G_real,
               G_im=G_im,lc_1_ref=lc_1_ref,lc_2_ref=lc_2_ref,data_2=data_2,n=n,m=m,fmin=fmin,fmax=fmax,seg_length=seg_length,bin_length=bin_length,GTI=GTI,coherence_corrector=False,spurious_sub=False)
    #print('mod_min_array',mod_min_array)
        partial_dG_calc_VEC = np.vectorize(partial_dG_calc,otypes=[object])
        result= partial_dG_calc_VEC(mod_min_array, mod_max_array,lc_subject)

    #print('mod_max_array',mod_max_array)

    
    #vectorize real/im calculator
    #partial_dG_calc_VEC = np.vectorize(partial_dG_calc,otypes=[object])
    #print('moving to g calc')
    # Get the real and imaginary parts of the power spectrum for the current modulation range
    #result= partial_dG_calc_VEC(mod_min_array, mod_max_array,lc_subject)#,data_1, lc_ref,GTI,bin_length, seg_length, fmin, fmax,spur_sub)
    # = zip(*result)
    #return result