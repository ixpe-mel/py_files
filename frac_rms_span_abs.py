#Calculates the fractional rms of a range of modulation angles
import numpy as np
import frac_rms_abs as frms
from functools import partial
def frac_rms_span(G_real_span,G_im_span,cr_span):
    frms_partial=partial(frms.frac_rms)
    frms_vec=np.vectorize(frms_partial)
    frms_span=frms_vec(G_real_span,G_im_span,cr_span)
    return frms_span
    