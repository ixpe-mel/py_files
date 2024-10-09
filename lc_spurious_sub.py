#Creates a lightcurve from spurious stokes and subtracts it from the original lightcurve
import numpy as np
from stingray import Lightcurve
def lc_spur_sub(data,lc_countrate,mod_bin_number,av_mod,bin_length,lc):
    """ Subtract spurious polarisation from a lightcurve. """
    q_spur=data['QSP'] # per event spurious stokes parameters for mod bin of interest
    u_spur=data['USP']

    Q_spur=np.sum(q_spur) # equivelent to big q
    U_spur=np.sum(u_spur) #equivalent to big u
    I=len(q_spur) # equivelent to I or N

    Q_spur_norm=np.sum(q_spur)/I #define normalised spur stokes 
    U_spur_norm=np.sum(u_spur)/I
        
        
    spur_sub=((lc_countrate/mod_bin_number)* ((Q_spur_norm*np.cos((2*av_mod)))+(U_spur_norm*np.sin(((2*av_mod))))))*(bin_length) #assumes qsm,usm const over time
    spur_sub_counts=[spur_sub]*len(lc.time)
    lc_spur=Lightcurve(lc.time,spur_sub_counts)                                
    lc_counts_subtracted=lc.counts-spur_sub_counts #subtracting spur lc
    lc=Lightcurve(lc.time,lc_counts_subtracted)
    return lc