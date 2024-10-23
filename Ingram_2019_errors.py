import numpy as np
import math
import lc_spurious_sub as lcss
import importlib
importlib.reload(lcss)
from stingray import Lightcurve, Crossspectrum, AveragedCrossspectrum, Powerspectrum
#Calculating Ingram 2019 Errorbars on rms and phase
def ingr_2019_errs_cc(mod_min,mod_max,mod_bin_number,real_G,im_G,lc_subject,cs_ref_real_mean,data_2,ps_2_mean,n,m,fmin,fmax,seg_length,bin_length,GTI,spurious_sub):
   
        # Creating Cross spectrum 'ith'
    #print(GTI)
    av_mod=(mod_min+mod_max)/2


    data_cut_2=data_2[(mod_min<=data_2['PHI']) & (data_2['PHI']<=mod_max)]
    lc2_subject=Lightcurve.make_lightcurve(data_cut_2['TIME'],dt=bin_length,gti=GTI)
    lc2_subject.apply_gtis()
    lc2_subject_countrate2=lc2_subject.meanrate
    if spurious_sub==True:
        print('Spurious sub on lc2')
        lc2_subject=lcss.lc_spur_sub(data_cut_2,lc2_subject_countrate2,mod_bin_number,av_mod,bin_length,lc2_subject)
    #lc2_subject_countrate2=lc2_subject.meanrate
    else:
        lc2_subject=lc2_subject 
    #print('lc subject type',type(lc_subject))
    ps_1_subject=Powerspectrum.from_lightcurve(lc_subject,seg_length,norm='frac')
    ps_1_subject_av=ps_1_subject.power[(fmin<=ps_1_subject.freq) & (ps_1_subject.freq<=fmax)].mean()
    #print('ps_1_subject_av',ps_1_subject_av)
    del ps_1_subject


    cs_ith=AveragedCrossspectrum.from_lightcurve(lc_subject,lc2_subject,seg_length,norm='frac')
    cs_ith_real_mean=cs_ith.power.real[(fmin<=cs_ith.freq) & (cs_ith.freq<=fmax)].mean()
    #print('cs_ith_real_mean',cs_ith_real_mean)
    #cs_ith_im_mean=cs_ith.power.imag[(fmin<=cs_ith.freq) & (cs_ith.freq<=fmax)].mean()
    del cs_ith
    #print('real G',type(real_G))
    #print('im G',type(im_G))
    modulus_G=np.sqrt((real_G**2)+(im_G**2))
    cs_ith_new=modulus_G**2/(cs_ref_real_mean)
    #print('cs_ith_new',cs_ith_new)
    # Power spectrum from file 2
    

    #print('modulus G',modulus_G)
    #Catching if the coherence is below 0 or above 1 and correcting....
    #print('modulus_G',modulus_G)
    #print('cs_ref_real_mean',cs_ref_real_mean)
    #print('cs_ith_real_mean',cs_ith_real_mean)
    coherence=modulus_G**2/(cs_ref_real_mean*cs_ith_real_mean)
    #print('coherence',coherence)

    corrected_coherence = np.clip(coherence, 0, 1)
    #print('ps_2_mean',ps_2_mean)
    #print('cs_ith_real_mean',cs_ith_real_mean)
    #print('new_coherence',new_coherence)
    #print(ps_1_subject_av)
    #print('new_coherence',new_coherence)
    #dG_2=

    dG_original=np.sqrt(      (ps_2_mean / (2*n*m)) * (ps_1_subject_av - corrected_coherence*(cs_ith_real_mean) )     )
    #print('dG_original',dG_original)
    #dG_new_csith=np.sqrt(      (ps_2_mean / (2*n*m)) * (ps_1_subject_av - corrected_coherence*(cs_ith_new) ) )
    #print('dG_new_csith',dG_new_csith)
    #dG_nogamma=np.sqrt(    (ps_2_mean/(2*n*m))   *   (ps_1_subject_av - (modulus_G**2/cs_ref_real_mean)    ))    
    #print('dG_nogamma',dG_nogamma)

    #dG_adam=np.sqrt( (  (ps_2_mean*ps_1_subject_av) - ( (ps_2_mean/cs_ref_real_mean) * modulus_G**2   )       )    / (2*n*m)      )
    #print(dG_adam)


    
    return dG_original


#Assuming no coherence correction, calculate Ingram 2019 errorbars dG from 

def ingr_2019_errs(real_G,im_G,lc_subject,cs_ref_real_mean,ps_2_mean,n,m,fmin,fmax,seg_length,bin_length,GTI,spurious_sub):
    modulus_G=np.sqrt((real_G**2)+(im_G**2))
    # Creating Powerspectrum
    #print('lc subject type',type(lc_subject))
    ps_1_subject=Powerspectrum.from_lightcurve(lc_subject,seg_length,norm='frac')
    ps_1_subject_av=ps_1_subject.power[(fmin<=ps_1_subject.freq) & (ps_1_subject.freq<=fmax)].mean()
    del ps_1_subject
        
    dG=np.sqrt(    (ps_2_mean/(2*n*m))   *   (ps_1_subject_av - (modulus_G**2/cs_ref_real_mean)    ))
    return dG

    #if square bracket negative set poisson noise 0 so its just ps1sub-/...