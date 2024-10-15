import numpy as np
from stingray import Lightcurve, Powerspectrum, AveragedCrossspectrum
import Ingram_2019_errors as I_19errs
import importlib
import lc_spurious_sub as lcss
importlib.reload(I_19errs)
def dG_calc(mod_min,mod_max,lc_subject,G_real,G_im,
            lc_1_ref,lc_2_ref,data_2,n,m,fmin,fmax,seg_length,bin_length,mod_bin_number,GTI,coherence_corrector=True,spurious_sub=True):
    #print(lc_subject)
    av_mod=(mod_min+mod_max)/2
    #lightcurve_2=Lightcurve.make_lightcurve(data_2['TIME'],dt=bin_length,gti=GTI)
    #lightcurve_2.apply_gtis()
    #lightcurve_2_countrate2=lightcurve_2.meanrate
    #lightcurve_2=lcss.lc_spur_sub(data_2,lightcurve_2_countrate2,mod_bin_number,av_mod,bin_length,lightcurve_2)
    
    
    ps_2=Powerspectrum.from_lightcurve(lc_2_ref,seg_length,norm='frac')
    ps_2_mean=ps_2.power[(fmin<=ps_2.freq) & (ps_2.freq<=fmax)].mean()
    print('ps_2_mean',ps_2_mean)

    cs_ref=AveragedCrossspectrum.from_lightcurve(lc_1_ref,lc_2_ref,seg_length,norm='frac')
    cs_ref_real_mean=cs_ref.power.real[(fmin<=cs_ref.freq) & (cs_ref.freq<=fmax)].mean()
    print('cs_ref_real_mean',cs_ref_real_mean)

    if coherence_corrector==True and spurious_sub==True:
        print('Applying coherence correction and spurious sub...')
        #print('Applying coherence correction...')
        dG=I_19errs.ingr_2019_errs_cc(mod_min,mod_max,mod_bin_number,G_real,G_im,lc_subject,cs_ref_real_mean,data_2,ps_2_mean,n,m,fmin,fmax,seg_length,bin_length,GTI,spurious_sub=True)
    elif coherence_corrector==True and spurious_sub==False:
        print('Applying coherence correction BUT NO SPUR SUB...')
        dG=I_19errs.ingr_2019_errs_cc(mod_min,mod_max,mod_bin_number,G_real,G_im,lc_subject,cs_ref_real_mean,data_2,ps_2_mean,n,m,fmin,fmax,seg_length,bin_length,GTI,spurious_sub=False)

    else:
        print('Applying no coherence correction or spur sub...')
        dG=I_19errs.ingr_2019_errs(real_G=G_real,im_G=G_im,lc_subject=lc_subject,cs_ref_real_mean=cs_ref_real_mean,ps_2_mean=ps_2_mean,n=n,m=m,fmin=fmin,fmax=fmax,seg_length=seg_length,bin_length=bin_length,GTI=GTI,spurious_sub=False)

    print('dG',dG)
        #dG=I_19errs.ingr_2019_errs(G_real,G_im,lc_subject,ps_2_mean,n,m,fmin,fmax,seg_length,bin_length,GTI)
    return dG
