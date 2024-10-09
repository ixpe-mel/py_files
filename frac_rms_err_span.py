import frac_rms_err as frms_err
import numpy as np
from functools import partial
def frac_rms_err_span(dG,norm_factor):
    frms_partial=partial(frms_err.frac_rms_err,norm_factor=norm_factor)
    frms_err_vec=np.vectorize(frms_partial)
    frms_err_span=frms_err_vec(dG)
    return frms_err_span