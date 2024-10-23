#New dG span that doesnt die when it reaches np.vectorize


import numpy as np
import dG_calc as dGc
from functools import partial
import importlib
importlib.reload(dGc)

def dG_span(mod_bin_number,
            G_real_span,G_im_span,
            lc_subject_span,lc_1_ref,lc_2_ref,
            data_2,
            n,m,fmin,fmax,seg_length,bin_length,GTI,ps_2_mean,cs_ref_real_mean,
            coherence_corrector=True,spurious_sub=True):

    #create modulation angle bins
    mod_min_global = np.radians(-90)
    mod_max_global = np.radians(90)
    aspace = np.linspace(mod_min_global, mod_max_global, mod_bin_number + 1)

    mod_min_array = aspace[:-1]
    mod_max_array = aspace[1:]
    av_mod = (mod_min_array + mod_max_array) / 2
    av_mod_err = (mod_max_array - mod_min_array) / 2
    result=[]
    if coherence_corrector==True and spurious_sub==True:
        #print('Applying coherence correction and spurious sub...')
        partial_dG_calc=partial(dGc.dG_calc
                                ,lc_1_ref=lc_1_ref,lc_2_ref=lc_2_ref,data_2=data_2,n=n,m=m,fmin=fmin,fmax=fmax,seg_length=seg_length,bin_length=bin_length,mod_bin_number=mod_bin_number,GTI=GTI,ps_2_mean=ps_2_mean,cs_ref_real_mean=cs_ref_real_mean,coherence_corrector=True,spurious_sub=True)

        
        results = [partial_dG_calc(min_val, max_val, span,real,im) for min_val, max_val, span ,real,im in zip(mod_min_array, mod_max_array, lc_subject_span,G_real_span,G_im_span)]
        #print(results)

        
        #for i in zip(mod_min_array,mod_max_array,lc_subject_span,G_real_span,G_im_span):
         #   result.append(partial_dG_calc(i[0],i[1],i[2],i[3],i[4]))

    elif coherence_corrector==True and spurious_sub==False:
        #print('Applying coherence correction BUT NO SPUR SUB...')
        partial_dG_calc=partial(dGc.dG_calc,
                                lc_1_ref=lc_1_ref,lc_2_ref=lc_2_ref,data_2=data_2,n=n,m=m,fmin=fmin,fmax=fmax,seg_length=seg_length,bin_length=bin_length,mod_bin_number=mod_bin_number,GTI=GTI,ps_2_mean=ps_2_mean,cs_ref_real_mean=cs_ref_real_mean,spurious_sub=False)

        #for i in zip(mod_min_array,mod_max_array,lc_subject_span,G_real_span,G_im_span):
        #    result.append(partial_dG_calc(i[0],i[1],i[2],i[3],i[4]))

        result = [partial_dG_calc(min_val, max_val, span,real,im) for min_val, max_val, span ,real,im in zip(mod_min_array, mod_max_array, lc_subject_span,G_real_span,G_im_span)]

    elif coherence_corrector==False and spurious_sub==True:
            
            partial_dG_calc=partial(dGc.dG_calc,
                                    lc_1_ref=lc_1_ref,lc_2_ref=lc_2_ref,data_2=data_2,n=n,m=m,fmin=fmin,fmax=fmax,seg_length=seg_length,bin_length=bin_length,mod_bin_number=mod_bin_number,GTI=GTI,ps_2_mean=ps_2_mean,cs_ref_real_mean=cs_ref_real_mean,coherence_corrector=False,spurious_sub=True)
    
            #for i in zip(mod_min_array,mod_max_array,lc_subject_span,G_real_span,G_im_span):
            #    result.append(partial_dG_calc(i[0],i[1],i[2],i[3],i[4]))

            result = [partial_dG_calc(min_val, max_val, span,real,im) for min_val, max_val, span ,real,im in zip(mod_min_array, mod_max_array, lc_subject_span,G_real_span,G_im_span)]

    else:
        partial_dG_calc=partial(dGc.dG_calc, 
                                lc_1_ref=lc_1_ref,lc_2_ref=lc_2_ref,data_2=data_2,n=n,m=m,fmin=fmin,fmax=fmax,seg_length=seg_length,bin_length=bin_length,mod_bin_number=mod_bin_number,GTI=GTI,ps_2_mean=ps_2_mean,cs_ref_real_mean=cs_ref_real_mean,coherence_corrector=False,spurious_sub=False)
        #for i in zip(mod_min_array,mod_max_array,lc_subject_span,G_real_span,G_im_span):
        #    result.append(partial_dG_calc(i[0],i[1],i[2],i[3],i[4]))
        result = [partial_dG_calc(min_val, max_val, span,real,im) for min_val, max_val, span ,real,im in zip(mod_min_array, mod_max_array, lc_subject_span,G_real_span,G_im_span)]