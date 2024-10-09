#Calculates error associated with phase in units of cycles
import numpy as np
def phase_err(G_real,G_im,dG):
    return dG / ( np.sqrt(G_real**2+G_im**2) * 2*np.pi )