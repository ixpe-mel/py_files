import numpy as np
import math
import importlib
from stingray import Powerspectrum

#Calculating Ingram 2019 Errorbars on rms and phase
def ingr_2019_errs_cc(real_G,im_G,lc_subject,cs_ref_real_mean,ps_2_mean,n,m,fmin,fmax,seg_length,stack,span2,n2,m2):
   
    print(lc_subject)
    ps_1_subject=Powerspectrum.from_lightcurve(lc_subject,seg_length,norm='abs')
    ps_1_subject_av=ps_1_subject.power[(fmin<=ps_1_subject.freq) & (ps_1_subject.freq<=fmax)].mean()

    #del ps_1_subject

    if stack==True:
        print(span2)
        ps_2_subject=Powerspectrum.from_lightcurve(span2,seg_length,norm='abs')
        ps_2_subject_av=ps_2_subject.power[(fmin<=ps_2_subject.freq) & (ps_2_subject.freq<=fmax)].mean()
        #del ps_2_subject

        ps_subject_stack=( (ps_1_subject.m*ps_1_subject_av) +(ps_2_subject.m*ps_2_subject_av) )/ (ps_1_subject.m + ps_2_subject.m)
        #cs_ref_real_mean=cs_ref_real_mean_obs2

        modulus_G=np.sqrt((real_G**2)+(im_G**2))

        cs_ith_new=modulus_G**2/(cs_ref_real_mean)

        coherence=modulus_G**2/(cs_ref_real_mean*cs_ith_new)

        corrected_coherence = np.clip(coherence, 0, 1)
    
        dG=np.sqrt(      (ps_2_mean / (2*(n+n2)*(m+m2))) * (ps_subject_stack - (corrected_coherence*cs_ith_new) )     )
    else:
        dG=np.sqrt(      (ps_2_mean / (2*n*m)) * (ps_1_subject_av - (corrected_coherence*cs_ith_new) )     )

    if math.isnan(dG):

        dG=np.sqrt(( (ps_2_mean*ps_1_subject_av) - modulus_G**2)/(2*n*m))
        #print('dG',dG)
    else:

        dG=dG


    
    return dG


#Assuming no coherence correction, calculate Ingram 2019 errorbars dG from 

def ingr_2019_errs(real_G,im_G,lc_subject,cs_ref_real_mean,ps_2_mean,n,m,fmin,fmax,seg_length,stack,span2,n2,m2):
    modulus_G=np.sqrt((real_G**2)+(im_G**2))
    # Creating Powerspectrum
    print(lc_subject)
    ps_1_subject=Powerspectrum.from_lightcurve(lc_subject,seg_length,norm='abs')
    ps_1_subject_av=ps_1_subject.power[(fmin<=ps_1_subject.freq) & (ps_1_subject.freq<=fmax)].mean()
    if stack==True:
        

        ps_2_subject=Powerspectrum.from_lightcurve(span2,seg_length,norm='abs')
        ps_2_subject_av=ps_2_subject.power[(fmin<=ps_2_subject.freq) & (ps_2_subject.freq<=fmax)].mean()
        #del ps_2_subject


        print(n)
        print(n2)
        print(m)
        print(m2)
        nstack=n+n2[0]
        mstack=m+m2[0]


        ps_subject_stack=np.array(( (ps_1_subject.m*ps_1_subject_av) +(ps_2_subject.m*ps_2_subject_av) )/ (ps_1_subject.m + ps_2_subject.m))
        dG=np.sqrt(    (np.array(ps_2_mean)/(2*(nstack*mstack)))   *   (ps_subject_stack - (np.array(modulus_G**2)/np.array(cs_ref_real_mean))    ))
    else:
        dG=np.sqrt(    (ps_2_mean/(2*n*m))   *   (ps_1_subject_av - (modulus_G**2/cs_ref_real_mean)    ))

    if math.isnan(dG):
        #print('dG is nan')
        dG=np.sqrt(((ps_2_mean*ps_1_subject_av) - modulus_G**2)/(2*n*m))
    else:
        dG=dG
    
    


    #print('dG=',dG)
    
    return dG