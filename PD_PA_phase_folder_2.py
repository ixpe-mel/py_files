import numpy as np
from astropy.io import fits
import math
import matplotlib.pyplot as plt
from stingray import Lightcurve
import fits_filter_gti as ffg


def PD_PA(file_1,file_2,file_3,ixpe_master_GTI,
           response_dir,Pmin,Pmax,phase_bin_number,v,v_dot,v_double_dot,t_0,output_file):
    
    GTI_master=list(np.loadtxt(str(ixpe_master_GTI)))

    #Loacing and cleaning data
    data_1,header_1,*_=lac.load_and_clean(file_1,Pmin,Pmax)
    data_2,header_2,*_=lac.load_and_clean(file_2,Pmin,Pmax)
    data_3,header_3,*_=lac.load_and_clean(file_3,Pmin,Pmax)


    data_1=ffg.filter_fits_gti(data_1,ixpe_master_GTI)
    data_2=ffg.filter_fits_gti(data_2,ixpe_master_GTI)
    data_3=ffg.filter_fits_gti(data_3,ixpe_master_GTI)

    #ignoring the response files for now till adam gets back to me

    