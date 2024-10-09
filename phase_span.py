#Vectorising phase caclulation over a set of calculations from different modulation angle values


import phase as ph
import numpy as np
def phase_span(G_real_span,G_im_span):
    phase_vec=np.vectorize(ph.phase)
    phase_span=phase_vec(G_real_span,G_im_span)
    return phase_span