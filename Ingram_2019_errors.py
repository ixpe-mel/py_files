import numpy as np
import math
import importlib
from stingray import Powerspectrum

#Calculating Ingram 2019 Errorbars on rms and phase
def ingr_2019_errs_cc(real_G,im_G,lc_subject,cs_ref_real_mean,ps_2_mean,n,m,fmin,fmax,seg_length):
   

    ps_1_subject=Powerspectrum.from_lightcurve(lc_subject,seg_length,norm='frac')
    ps_1_subject_av=ps_1_subject.power[(fmin<=ps_1_subject.freq) & (ps_1_subject.freq<=fmax)].mean()

    del ps_1_subject



    modulus_G=np.sqrt((real_G**2)+(im_G**2))

    cs_ith_new=modulus_G**2/(cs_ref_real_mean)

    coherence=modulus_G**2/(cs_ref_real_mean*cs_ith_new)

    corrected_coherence = np.clip(coherence, 0, 1)
  
    dG=np.sqrt(      (ps_2_mean / (2*n*m)) * (ps_1_subject_av - (corrected_coherence*cs_ith_new) )     )
    print('Pr3=',ps_2_mean)
    print('Ps12=',ps_1_subject_av)
    print('mod G=',modulus_G)
#    print('Pco_r_123=',cs_ref_real_mean)

    #print('dG_original',dG_original)
    #dG=np.sqrt(      (ps_2_mean / (ps_1_subject_av - corrected_coherence*(cs_ith_new) ) )/(2*n*m)  )
    #dG=np.sqrt( ((ps_2_mean*ps_1_subject_av)-((ps_2_mean*modulus_G**2)/cs_ref_real_mean))  /  (2*n*m) )
    #print('new dg',dG)
     #print('dG_new_csith',dG_new_csith)
    #dG_nogamma=np.sqrt(    (ps_2_mean/(2*n*m))   *   (ps_1_subject_av - (modulus_G**2/cs_ref_real_mean)    ))    
#    dG=np.sqrt( (  (ps_2_mean*ps_1_subject_av) - ( (ps_2_mean/cs_ref_real_mean) * modulus_G**2   )       )    / (2*n*m)      )

    if math.isnan(dG):
#        print('Poisson noise is now 0')
        #print('ps_2_mean',ps_2_mean)
        #print('ps_1_subject_av',ps_1_subject_av)
 #       print('modulus_G',modulus_G)
 #       print('1st bit',(ps_2_mean*ps_1_subject_av))
 #       print('2nd bit',(modulus_G**2))
        dG=np.sqrt(( (ps_2_mean*ps_1_subject_av) - modulus_G**2)/(2*n*m))
 #       print('dG',dG)
    else:
#        print('Catch not used')
#        print('Pco_r_123=',cs_ref_real_mean)
        dG=dG
    print('dG',dG)
    #print(dG_adam)


    
    return dG


#Assuming no coherence correction, calculate Ingram 2019 errorbars dG from 

def ingr_2019_errs(real_G,im_G,lc_subject,cs_ref_real_mean,ps_2_mean,n,m,fmin,fmax,seg_length):
    modulus_G=np.sqrt((real_G**2)+(im_G**2))
    # Creating Powerspectrum

    ps_1_subject=Powerspectrum.from_lightcurve(lc_subject,seg_length,norm='frac')
    ps_1_subject_av=ps_1_subject.power[(fmin<=ps_1_subject.freq) & (ps_1_subject.freq<=fmax)].mean()
    del ps_1_subject
        
    dG=np.sqrt(    (ps_2_mean/(2*n*m))   *   (ps_1_subject_av - (modulus_G**2/cs_ref_real_mean)    ))
    print('Pr3=',ps_2_mean)
    print('Ps12=',ps_1_subject_av)
    print('mod G=',modulus_G)
    print('Pco_r_123=',cs_ref_real_mean)
    if math.isnan(dG):
#        print('dG is nan')
        dG=np.sqrt(((ps_2_mean*ps_1_subject_av) - modulus_G**2)/(2*n*m))
    else:
        dG=dG
    
    
    print('dG=',dG)
    
    return dG
