#Calculating Stokes parameters using the IXPEOBSSIM convention

import numpy as np
def calculate_stokes(W_MOM_phase_cut, q_phase_cut, u_phase_cut, Aeff_event, Aeff_mu_event):
    

    #print('inside wmom type',type(W_MOM_phase_cut))
    #print('inside aeff type',Aeff_event)

    I_NEFF = np.sum(W_MOM_phase_cut/Aeff_event)
    #print('I_NEFF',I_NEFF)
    dI_NEFF_sqrd = np.sum((W_MOM_phase_cut / Aeff_event) ** 2)

    Q_NEFF = np.sum((W_MOM_phase_cut * q_phase_cut ) / Aeff_mu_event)
    dQ_NEFF_sqrd = np.sum(((W_MOM_phase_cut * q_phase_cut) / Aeff_event) ** 2)

    U_NEFF = np.sum((W_MOM_phase_cut * u_phase_cut ) / Aeff_mu_event)
    dU_NEFF_sqrd = np.sum(((W_MOM_phase_cut * u_phase_cut) / Aeff_event) ** 2)

    return I_NEFF, np.sqrt(dI_NEFF_sqrd), Q_NEFF, np.sqrt(dQ_NEFF_sqrd), U_NEFF, np.sqrt(dU_NEFF_sqrd)