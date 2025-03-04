#Modelling the cross spectrum t be used in measuring Q,U over nu
import numpy as np

def cross_spec_model_real(J,phi,A,B,C):
    Re_G=(1/J) * ( A + (B*np.cos(2*phi)) + (C*np.sin(2*phi)) )
    return Re_G

#The imaginary sinusoid does not have the A term
def cross_spec_model_imag(J,phi,B,C):
    Im_G=(1/J) * ( (B*np.sin(2*phi)) + (C*np.cos(2*phi)) )
    return Im_G


#Null hypothesis model (ie if PD/PA are constant then the cross spectrum becomes:)
def cross_spec_model_null(J,Q_norm,U_norm,phi,C_nu_mag_sqrd):
    ReG_null=(1/J) *C_nu_mag_sqrd* (1 + Q_norm*np.cos(2*phi) + U_norm*np.sin(2*phi))
    
    return ReG_null
