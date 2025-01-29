import numpy as np
from stingray import Lightcurve, AveragedCrossspectrum

#calculate the normalisation of fractional frac_rms
def rms_normalisation(data_12,data_3,bin_length,seg_length,fmin,fmax,gti,norm='frac'):
        
        GTI=list(np.loadtxt(str(gti)))
        data_12_times = data_12['TIME']
        
        lightcurve_12=Lightcurve.make_lightcurve(data_12_times,dt=bin_length,gti=GTI)
        lightcurve_12.apply_gtis()

        lc_ref=Lightcurve.make_lightcurve(data_3['TIME'],dt=bin_length,gti=GTI)
        lc_ref.apply_gtis()

        avg_cs = AveragedCrossspectrum.from_lightcurve(lightcurve_12,lc_ref,seg_length,norm=norm)
        del lightcurve_12

        av_power_norm_real=avg_cs.power.real[(fmin<=avg_cs.freq) & (avg_cs.freq<=fmax)].mean()
        norm_factor=(np.sqrt((fmax-fmin))/np.sqrt(av_power_norm_real))
        del avg_cs

        return norm_factor,lc_ref,GTI

