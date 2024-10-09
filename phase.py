
import numpy as np
def phase(G_real,G_im):
    return np.arctan(G_im/G_real)/(2*np.pi)

def phase_err(G_real,G_im,dG):
    return dG/(abs(G_real,G_im)*2*np.pi)