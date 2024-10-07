import numpy as np
import math

from stingray import Lightcurve, Crossspectrum, AveragedCrossspectrum, Powerspectrum
#Calculating Ingram 2019 Errorbars on rms and phase
def ingr_2019_errs_cc(mod_min,mod_max,real_G,im_G,lc_subject,data_2,ps_2_mean,n,m,fmin,fmax,seg_length,bin_length,GTI):
        # Creating Cross spectrum 'ith'

    
    data_cut_2=data_2[(mod_min<=data_2['PHI']) & (data_2['PHI']<=mod_max)]
    lc2_subject=Lightcurve.make_lightcurve(data_cut_2['TIME'],dt=bin_length,gti=GTI)
    lc2_subject.apply_gtis()
    lc2_subject_countrate2=lc2_subject.meanrate

    ps_1_subject=Powerspectrum.from_lightcurve(lc_subject,seg_length,norm='frac')
    ps_1_subject_av=ps_1_subject.power[(fmin<=ps_1_subject.freq) & (ps_1_subject.freq<=fmax)].mean()
    del ps_1_subject


    cs_ith=AveragedCrossspectrum.from_lightcurve(lc_subject,lc2_subject,seg_length,norm='frac')
    cs_ith_real_mean=cs_ith.power.real[(fmin<=cs_ith.freq) & (cs_ith.freq<=fmax)].mean()
    cs_ith_im_mean=cs_ith.power.imag[(fmin<=cs_ith.freq) & (cs_ith.freq<=fmax)].mean()
    del cs_ith
    # Power spectrum from file 2
    modulus_G=np.sqrt((real_G**2)+(im_G**2))


    #Catching if the coherence is below 0 or above 1 and correcting....

    coherence=np.array(modulus_G**2/real_G*cs_ith_real_mean)
    print('coherence',coherence)

    new_coherence = np.clip(coherence, 0, 1)
    print('new_coherence',new_coherence)
    dG=np.sqrt(      (ps_2_mean / (n*m)) * (ps_1_subject_av - new_coherence*(cs_ith_real_mean) )     )

    return dG


#Assuming no coherence correction, calculate Ingram 2019 errorbars dG from 

def ingr_2019_errs(real_G,im_G,lc_subject,ps_2_mean,n,m,fmin,fmax,seg_length,bin_length,GTI):
    modulus_G=np.sqrt((real_G**2)+(im_G**2))
    # Creating Powerspectrum
        
    ps_1_subject=Powerspectrum.from_lightcurve(lc_subject,seg_length,norm='frac')
    ps_1_subject_av=ps_1_subject.power[(fmin<=ps_1_subject.freq) & (ps_1_subject.freq<=fmax)].mean()
    del ps_1_subject
        
    dG=np.sqrt((1/(2*(n*m)))*((ps_1_subject_av*ps_2_mean)-((ps_2_mean/(real_G))*(modulus_G**2))))
    return dG