#Calculates the fractional rms of a range of modulation angles
import numpy as np
import frac_rms as frms
from functools import partial
def frac_rms_span(G_real_span,G_im_span,norm_factor):
    frms_partial=partial(frms.frac_rms,norm_factor=norm_factor)
    frms_vec=np.vectorize(frms_partial)
    frms_span=frms_vec(G_real_span,G_im_span)
    return frms_span
    