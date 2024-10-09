#Calculates the phase error from a span of G,dG

import phase_err as p_err
import numpy as np
from functools import partial
def phase_err_span(G_real,G_im,dG):
    phase_err_vec=np.vectorize(p_err.phase_err)
    phase_err_span=phase_err_vec(G_real,G_im,dG)
    return phase_err_span