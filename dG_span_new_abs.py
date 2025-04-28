
import numpy as np
import dG_calc_abs as dGc
from functools import partial
import importlib
importlib.reload(dGc)


def dG_span(G_real_span,G_im_span,
            lc_subject_span,n_span,m_span,
            fmin,fmax,seg_length,ps_2_mean,cs_ref_real_mean,
            coherence_corrector,stack,span2=None,n2=None,m2=None):


 #   result=[]
 #   if coherence_corrector==True:
 #       print(dGc.dG_calc.__code__.co_varnames)

#        partial_dG_calc=partial(dGc.dG_calc,n=n_span[0],m=m_span[0],fmin=fmin,fmax=fmax,
#                                seg_length=seg_length,
#                                ps_2_mean=ps_2_mean,cs_ref_real_mean=cs_ref_real_mean,coherence_corrector=True,stack=stack,lc_subject_obs2=lc_subject_obs2,n2=n2,m2=m2)
#        result = [partial_dG_calc(span,real,im,span2) for span ,real,im, span2  in zip( lc_subject_span,G_real_span,G_im_span,lc_subject_obs2)]

#    else:
#        print(dGc.dG_calc.__code__.co_varnames)

#        partial_dG_calc=partial(dGc.dG_calc, n=n_span[0],m=m_span[0],fmin=fmin,fmax=fmax,
#                                seg_length=seg_length,ps_2_mean=ps_2_mean,
#                                cs_ref_real_mean=cs_ref_real_mean,
#                                coherence_corrector=False,stack=stack,lc_subject_obs2=lc_subject_obs2,n2=n2,m2=m2)
#        
#        result = [partial_dG_calc( span,real,im,span2) for span ,real,im,span2 in zip( lc_subject_span,G_real_span,G_im_span,lc_subject_obs2)]
#    return result

    if coherence_corrector == True:
        print(dGc.dG_calc.__code__.co_varnames)

        def partial_dG_calc(span, real, im, span2):
            return dGc.dG_calc(span, real, im,
                            n=n_span[0], m=m_span[0],
                            fmin=fmin, fmax=fmax,
                            seg_length=seg_length,
                            ps_2_mean=ps_2_mean,
                            cs_ref_real_mean=cs_ref_real_mean,
                            coherence_corrector=True, stack=stack,span2=span2,
                            n2=n2, m2=m2)

    else:
        print(dGc.dG_calc.__code__.co_varnames)

        def partial_dG_calc(span, real, im, span2):
            return dGc.dG_calc(span, real, im,
                            n=n_span[0], m=m_span[0],
                            fmin=fmin, fmax=fmax,
                            seg_length=seg_length,
                            ps_2_mean=ps_2_mean,
                            cs_ref_real_mean=cs_ref_real_mean,
                            coherence_corrector=False, stack=stack,span2=span2,
 
                            n2=n2, m2=m2)

    result = [partial_dG_calc(span, real, im, span2) for span, real, im, span2 in zip(lc_subject_span, G_real_span, G_im_span, span2)]
    return result
#    result = [partial_dG_calc(span, real, im, span2) for span ,real,im, span2  in zip( lc_subject_span,G_real_span,G_im_span,lc_subject_obs2)]