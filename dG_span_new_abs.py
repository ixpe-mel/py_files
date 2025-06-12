
import numpy as np
import dG_calc_abs as dGc
from functools import partial
import importlib
importlib.reload(dGc)


def dG_span(G_real_span,G_im_span,
            lc_subject_span,n_span,m_span,
            fmin,fmax,seg_length,ps_2_mean,cs_ref_real_mean,
            coherence_corrector):


    if coherence_corrector == True:
        #print(dGc.dG_calc.__code__.co_varnames)

        def partial_dG_calc(span, real, im):
            return dGc.dG_calc(span, real, im,
                            n=n_span[0], m=m_span[0],
                            fmin=fmin, fmax=fmax,
                            seg_length=seg_length,
                            ps_2_mean=ps_2_mean,
                            cs_ref_real_mean=cs_ref_real_mean,
                            coherence_corrector=True)

    else:
        #print(dGc.dG_calc.__code__.co_varnames)

        def partial_dG_calc(span, real, im):
            return dGc.dG_calc(span, real, im,
                            n=n_span[0], m=m_span[0],
                            fmin=fmin, fmax=fmax,
                            seg_length=seg_length,
                            ps_2_mean=ps_2_mean,
                            cs_ref_real_mean=cs_ref_real_mean,
                            coherence_corrector=False)
    
    result = [partial_dG_calc(span, real, im) for span, real, im in zip(lc_subject_span, G_real_span, G_im_span)]
    return result
#    result = [partial_dG_calc(span, real, im, span2) for span ,real,im, span2  in zip( lc_subject_span,G_real_span,G_im_span,lc_subject_obs2)]