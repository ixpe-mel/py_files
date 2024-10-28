#Generates a modulation angle based on a given modulation factor , PD and PA
#Parameters:
#CDF0: Initial CDF value ie the random number generated
#PD: PD value
#PA: PA value
#mu: Modulation factor
import numpy as np
from scipy import optimize
def modulation_angle_generator(CDF0,muPD,PA):
    phimax=np.pi/2
    
    phimin=-np.pi/2
    mod_function=(lambda phi: (1/np.pi) * (1 + (muPD * np.cos(2 * (PA-phi)) ) ))
    mod_function_derivative=(lambda phi: ((1/np.pi)*2*muPD*np.sin(PA-phi)))
    modmin=np.radians(-90)        
    initial_guess=(CDF0*np.pi)+modmin
    cdf_NR=lambda phi: ( ( (1/(2*np.pi)) * ( (2*phi) + np.pi + muPD * (-np.sin(2*PA) + np.sin( (2*phi)-(2*PA) ) ) ) ) - CDF0 )
    cdf_NR_PA0=lambda phi:( ( (1/(2*np.pi)) * ( (2*phi) + np.pi + muPD * ( + np.sin( (2*phi) ) ) ) ) - CDF0 )
    NR_root=optimize.newton(cdf_NR,x0=initial_guess,fprime=mod_function,fprime2=mod_function_derivative,tol=10e-4,maxiter=10000000000000000000)
    NR_root_PA0=optimize.newton(cdf_NR_PA0,x0=initial_guess,fprime=mod_function,fprime2=mod_function_derivative,tol=10e-10,maxiter=10000000000000000000)
    #NR_root= 2 * PA - NR_root
    
    NR_root_PA0= NR_root_PA0 + PA
    range_ = phimax - phimin
    while NR_root_PA0 > phimax:
        NR_root_PA0 -= range_
    while NR_root_PA0 < phimin:
        NR_root_PA0 += range_
    #print('Modulation angle calculated (deg):',np.degrees(NR_root_PA0))
    return NR_root_PA0
