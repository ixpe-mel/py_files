import numpy as np
import math
import importlib
from stingray import Powerspectrum

#Calculating Ingram 2019 Errorbars on rms and phase
def ingr_2019_errs_cc(real_G,im_G,lc_subject,cs_ref_real_mean,ps_2_mean,n,m,fmin,fmax,seg_length):
   
    print(lc_subject)
    ps_1_subject=Powerspectrum.from_lightcurve(lc_subject,seg_length,norm='abs')
    ps_1_subject_av=ps_1_subject.power[(fmin<=ps_1_subject.freq) & (ps_1_subject.freq<=fmax)].mean()

    #del ps_1_subject


    dG=np.sqrt(      (ps_2_mean / (2*n*m)) * (ps_1_subject_av - (corrected_coherence*cs_ith_new) )     )

    if math.isnan(dG):

        dG=np.sqrt(( (ps_2_mean*ps_1_subject_av) - modulus_G**2)/(2*n*m))
        #print('dG',dG)
    else:

        dG=dG


    
    return dG


#Assuming no coherence correction, calculate Ingram 2019 errorbars dG from 

def ingr_2019_errs(real_G,im_G,lc_subject,cs_ref_real_mean,ps_2_mean,n,m,fmin,fmax,seg_length):
    modulus_G=np.sqrt((real_G**2)+(im_G**2))
    # Creating Powerspectrum
    #print(lc_subject)
    ps_1_subject=Powerspectrum.from_lightcurve(lc_subject,seg_length,norm='abs')
    ps_1_subject_av=ps_1_subject.power[(fmin<=ps_1_subject.freq) & (ps_1_subject.freq<=fmax)].mean()

    dG=np.sqrt(    (ps_2_mean/(2*n*m))   *   (ps_1_subject_av - (modulus_G**2/cs_ref_real_mean)    ))

    if math.isnan(dG):
        #print('dG is nan')
        dG=np.sqrt(((ps_2_mean*ps_1_subject_av) - modulus_G**2)/(2*n*m))
    else:
        dG=dG
    
    


    #print('dG=',dG)
    
    return dG