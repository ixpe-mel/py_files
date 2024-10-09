#Calculates the error on frac rms from dG and normalisation factor

import numpy as np

def frac_rms_err(dG, norm_factor):
    frac_rms_err = dG*norm_factor
    return frac_rms_err