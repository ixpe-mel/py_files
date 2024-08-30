#Calculating the phase of a photon on its spin period
import numpy as np
import math

def phase(v,t,t_0):
    phi=(v)*(t-t_0)
    #print(phi)
    floor_vec=np.vectorize(math.floor)
    phi=phi-floor_vec(phi)
    #phi=2*np.pi*phi
    #print(phi)
    return phi
    
def phase_v_dot(v,v_dot,t,t_0):
    #print(v_dot)
    phi=( v*(t-t_0) ) + (np.pi*v_dot*(t-t_0)**2)
    phi=phi%1
    return phi
    
def phase_v_double_dot(v,v_dot,v_double_dot,t,t_0):
    #print(v_dot)
    phi=( (v*(t-t_0)) + (1/2)*(v_dot)*(t-t_0)**2 + (1/6)*(v_double_dot)*(t-t_0)**3 )
    phi=phi%1
    return phi    
 #v is the taylor expansion of v(t) approx v|t=to +vdot|t=to *(t-t0) +vdoubledot/2 * (t-t0)**2
 
 #so phase becomes
 
 # phase(t)=2pi{v*(t-to)+vdot/2 * (t-to)**2 + vdoubledot/6 * (t-to)**3}
   
