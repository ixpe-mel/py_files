import numpy as np
from stingray import Powerspectrum
import Ingram_2019_errors as I_19errs
def dG_calc(mod_min,mod_max,G_real,G_im,lc_subject,
            lc_ref,data_2,n,m,fmin,fmax,seg_length,bin_length,GTI,coherence_corrector=True):

    ps_2=Powerspectrum.from_lightcurve(lc_ref,seg_length,norm='frac')
    ps_2_mean=ps_2.power[(fmin<=ps_2.freq) & (ps_2.freq<=fmax)].mean()

    if coherence_corrector==True:
        print('Applying coherence correction...')
        dG=I_19errs.ingr_2019_errs_cc(mod_min,mod_max,G_real,G_im,lc_subject,data_2,ps_2_mean,n,m,fmin,fmax,seg_length,bin_length,GTI)
    else:
        print('Applying no coherence correction...')
        dG=I_19errs.ingr_2019_errs(real_G=G_real,im_G=G_im,lc_subject=lc_subject,ps_2_mean=ps_2_mean,n=n,m=m,fmin=fmin,fmax=fmax,seg_length=seg_length,bin_length=bin_length,GTI=GTI)

        #dG=I_19errs.ingr_2019_errs(G_real,G_im,lc_subject,ps_2_mean,n,m,fmin,fmax,seg_length,bin_length,GTI)
    return dG
