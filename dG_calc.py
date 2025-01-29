import numpy as np
from stingray import Lightcurve, Powerspectrum, AveragedCrossspectrum
import Ingram_2019_errors as I_19errs
import importlib
import lc_spurious_sub as lcss
importlib.reload(I_19errs)
def dG_calc(lc_subject,G_real,G_im,n,m,fmin,fmax,seg_length,ps_2_mean,cs_ref_real_mean,coherence_corrector,norm):

    if coherence_corrector==True:
        #print('Applying coherence correction')
        dG=I_19errs.ingr_2019_errs_cc(G_real,G_im,lc_subject,cs_ref_real_mean,ps_2_mean,n,m,fmin,fmax,seg_length,norm)

    else:
        #print('Applying no coherence correction...')
        dG=I_19errs.ingr_2019_errs(real_G=G_real,im_G=G_im,lc_subject=lc_subject,cs_ref_real_mean=cs_ref_real_mean,ps_2_mean=ps_2_mean,n=n,m=m,fmin=fmin,fmax=fmax,seg_length=seg_length,norm=norm)

    return dG
