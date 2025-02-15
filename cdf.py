#Generates a modulation angle based on a given modulation factor , PD and PA

#Parameters:
#CDF0: Initial CDF value ie the random number generated
#PD: PD value
#PA: PA value
#mu: Modulation factor

import numpy as np
from scipy import optimize



def modulation_angle_generator(CDF0,muPD,PA,phimin,phimax):#,spur_sub,qsp,usp,mue):
    #loop_count=
    #print('inside modulation_angle_generator')


    
    mod_function=(lambda phi: (1/np.pi) * (1 + (muPD * np.cos(2 * (PA-phi)) ) ))
    mod_function_derivative=(lambda phi: ((1/np.pi)*2*muPD*np.sin(PA-phi)))
    #print('mod function assinged')
    #modmin=np.radians(-90)        
    initial_guess=(CDF0*np.pi)+phimin
    #cdf_NR=lambda phi: ( ( (1/(2*np.pi)) * ( (2*phi) + np.pi + muPD * (-np.sin(2*PA) + np.sin( (2*phi)-(2*PA) ) ) ) ) - CDF0 )
    cdf_NR_PA0=lambda phi:( ( (1/(2*np.pi)) * ( (2*phi) + np.pi + muPD * ( + np.sin( (2*phi) ) ) ) ) - CDF0 )
    #NR_root=optimize.newton(cdf_NR,x0=initial_guess,fprime=mod_function,fprime2=mod_function_derivative,tol=10e-4,maxiter=10000000000000000000)
    NR_root_PA0,rr=optimize.newton(cdf_NR_PA0,x0=initial_guess,fprime=mod_function,fprime2=mod_function_derivative,tol=10e-6,maxiter=100000,disp=False,full_output=True)
    #print(rr.converged)
    #NR_root_PA0=optimize.newton(cdf_NR_PA0,x0=initial_guess,fprime=mod_function,fprime2=mod_function_derivative,tol=10e-1,maxiter=10000000)
    #print(NR_root_PA0)
    #print(result_conv)
    #NR_root_PA0=NR_root_PA0[0]
    #NR_root_conv=NR_root_PA0[1]
    #print('root found')
    #NR_root= 2 * PA - NR_root
    #print('CDF0:',CDF0)
    NR_root_PA0= NR_root_PA0 + PA

    range_ = phimax - phimin
    while NR_root_PA0 > phimax:
        NR_root_PA0 -= range_
    while NR_root_PA0 < phimin:
        NR_root_PA0 += range_
    #print('muPD:',muPD)
    #print('Modulation angle calculated (deg):',np.degrees(NR_root_PA0))


    #mod_function=(lambda phi: (1/np.pi) * (1 + (muPD * np.cos(2 * (PA-phi)) ) ))
    #mod_function_derivative=(lambda phi: ((1/np.pi)*2*muPD*np.sin(PA-phi)))
    #print('mod function assinged')
    #modmin=np.radians(-90)        
    #initial_guess=(CDF0*np.pi)+phimin
    #cdf_NR=lambda phi: ( ( (1/(2*np.pi)) * ( (2*phi) + np.pi + muPD * (-np.sin(2*PA) + np.sin( (2*phi)-(2*PA) ) ) ) ) - CDF0 )
    #cdf_NR_PA0=lambda phi:( ( (1/(2*np.pi)) * ( (2*phi) + np.pi + muPD * ( + np.sin( (2*phi) ) ) ) ) - CDF0 )
    #NR_root=optimize.newton(cdf_NR,x0=initial_guess,fprime=mod_function,fprime2=mod_function_derivative,tol=10e-4,maxiter=10000000000000000000)
    #NR_root_PA0=optimize.newton(cdf_NR_PA0,x0=initial_guess,fprime=mod_function,fprime2=mod_function_derivative,tol=10e-10,maxiter=10000000000000000000)
    #NR_root_PA0=optimize.newton(cdf_NR_PA0,x0=initial_guess,fprime=mod_function,fprime2=mod_function_derivative,tol=10e-1,maxiter=10000000)
    #print(NR_root_PA0)
    #print(result_conv)
    #NR_root_PA0=NR_root_PA0[0]
    #NR_root_conv=NR_root_PA0[1]
    #print('root found')
    #NR_root= 2 * PA - NR_root
    #print('CDF0:',CDF0)
    #NR_root_PA0= NR_root_PA0 + PA

    #range_ = phimax - phimin
    #while NR_root_PA0 > phimax:
    #    NR_root_PA0 -= range_
    #while NR_root_PA0 < phimin:
    #    NR_root_PA0 += range_
#print(loop_count)
#print(NR_root_PA0)
#if loop_count % 1000 == 0:
#    print(f"Progress report: Completed {loop_count} iterations")
#print('Modulation angle calculated (deg):',np.degrees(NR_root_PA0))
#print(NR_root_PA0)
#np.save('/home/c2032014/converge_info.npy',NR_root_conv)
    return NR_root_PA0,rr

