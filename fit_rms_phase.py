#Takes rms and phase with corresponding errors and resturns the statistical significance of each fit
import sys
sys.path.append('/home/c2032014/py_files/')
import fit_models as fm
import chi_square as cs
import numpy as np
import scipy
import scipy.special
from scipy.optimize import curve_fit

# Fitting a straight line to rms/phase and calculating chi squared values
def fit_line(av_mod,frac_rms,frac_rms_err,phase,phase_err):    
    
    #Fitting line to rms
    parameters_line_fracrms, covariance_fracrms = curve_fit(fm.line, av_mod ,frac_rms,maxfev=50000) #defining the returns of curvefit
    fit_y_line_fracrms = fm.line(np.array(av_mod),fit_const_fracrms)
    #print('frac rms const fit {}'.format(fit_const_fracrms))
    
    
    #Fitting line to phase
    parameters_line_phase, covariance_phase = curve_fit(fm.line, av_mod ,phase,maxfev=50000) #defining the returns of curvefit
    fit_y_line_phase = fm.line(np.array(av_mod),fit_const_phase)
    
    dof_line_rms,dof_line_phase=len(av_mod)-1


    frac_rms_line_chi=cs.chi_square(frac_rms,fit_y_line_fracrms,frac_rms_err)#phase,phase_err,fit_y_line_phase,dof_line)
    phase_line_chi=cs.chi_square(phase,phase_err,fit_y_line_phase,dof_line_rms,dof_line_phase)

    reduced_chi_line=cs.reduced_chi_square(frac_rms_line_chi,phase_line_chi,dof_line)
    return parameters_line_fracrms,fit_y_line_fracrms,parameters_line_phase,fit_y_line_phase,dof_line,frac_rms_line_chi,phase_line_chi,reduced_chi_line

def fit_sine_90(av_mod,frac_rms,frac_rms_err,phase,phase_err):
    
    
    # Fitting 90 degree sine wave to rms/phase
    parameters_sin90_fracrms, sin90_covariance_fracrms = curve_fit(fm.sin_90, av_mod, frac_rms, maxfev=1000)#,p0=real_sin_initial_guess)
    parameters_sin90_phase, sin90_covariance_phase = curve_fit(fm.sin_90, av_mod, phase,maxfev=50000)#,p0=im_sin_initial_guess)


    fit_y_sin90_rms = fm.sin_90(av_mod, parameters_sin90_fracrms)#fit_A_0_90_rms, fit_A_90_rms,fit_delta_90_rms)
    fit_y_sin90_phase=fm.sin_90(av_mod,parameters_sin90_phase)#fit_A_0_90_phase,fit_A_90_phase,fit_delta_90_phase)


    #Calculating chi sqr value of 90 sine fit
    dof_sin90_rms=len(av_mod)-3
    dof_sin90_phase=len(av_mod)-3


    frac_rms_sin90_chi=cs.chi_square(frac_rms,fit_y_sin90_fracrms,frac_rms_err)#phase,phase_err,fit_y_line_phase,dof_line)
    phase_sin90_chi=cs.chi_square(phase,phase_err,fit_y_sin90_phase,dof_line)
    reduced_chi_sin90=cs.reduced_chi_square(frac_rms_sin90_chi,phase_sin90_chi,dof_line)
    return parameters_sin90_fracrms,fit_y_sin90_fracrms,parameters_sin90_phase,fit_y_sin90_phase,dof_sin90_rms,dof_sin90_phase,frac_rms_sin90_chi,phase_sin90_chi,reduced_chi_sin90

def fit_sine_180(av_mod,frac_rms,frac_rms_err,phase,phase_err):
 
    
    
    # Fitting 90 degree sine wave to rms/phase
    parameters_sin180_fracrms, sin180_covariance_fracrms = curve_fit(fm.sin_180, av_mod, frac_rms, maxfev=1000)#,p0=real_sin_initial_guess)
    parameters_sin180_phase, sin180_covariance_phase = curve_fit(fm.sin_180, av_mod, phase,maxfev=50000)#,p0=im_sin_initial_guess)


    fit_y_sin180_rms = fm.sin_180(av_mod, fit_A_0_180_rms, fit_A_180_rms,fit_delta_180_rms)
    fit_y_sin180_phase=fm.sin_180(av_mod,fit_A_0_180_phase,fit_A_180_phase,fit_delta_180_phase)


    #Calculating chi sqr value of 90 sine fit
    dof_sin180_rms=len(av_mod)-3
    dof_sin180_phase=len(av_mod)-3


    frac_rms_sin180_chi=cs.chi_square(frac_rms,fit_y_sin180_fracrms,frac_rms_err)#phase,phase_err,fit_y_line_phase,dof_line)
    phase_sin180_chi=cs.chi_square(phase,phase_err,fit_y_sin180_phase,dof_line)
    reduced_chi_sin180=cs.reduced_chi_square(frac_rms_sin180_chi,phase_sin180_chi,dof_sin180_rms,dof_sin180_phase)
    return parameters_sin180_fracrms,fit_y_sin180_fracrms,parameters_sin180_phase,fit_y_sin180_phase,dof_sin180_rms,dof_sin180_phase,frac_rms_sin180_chi,phase_sin180_chi,reduced_chi_sin180


def fit_sinsum(av_mod,frac_rms,frac_rms_err,phase,phase_err):
 
    
    
    # Fitting 90 degree sine wave to rms/phase
    parameters_sinsum_fracrms, sinsum_covariance_fracrms = curve_fit(fm.sinsum, av_mod, frac_rms, maxfev=1000)#,p0=real_sin_initial_guess)
    parameters_sinsum_phase, sinsum_covariance_phase = curve_fit(fm.sinsum, av_mod, phase,maxfev=50000)#,p0=im_sin_initial_guess)


    fit_y_sinsum_rms = fm.sinsum(av_mod,parameters_sinsum_fracrms)# fit_A_0_180_rms, fit_A_180_rms,fit_delta_180_rms)
    fit_y_sinsum_phase=fm.sin_180(av_mod,parameters_sinsum_phase)#fit_A_0_180_phase,fit_A_180_phase,fit_delta_180_phase)


    #Calculating chi sqr value of 90 sine fit
    dof_sinsum_rms=len(av_mod)-3
    dof_sinsum_phase=len(av_mod)-3


    frac_rms_sinsum_chi=cs.chi_square(frac_rms,fit_y_sinsum_fracrms,frac_rms_err)#phase,phase_err,fit_y_line_phase,dof_line)
    phase_sinsum_chi=cs.chi_square(phase,phase_err,fit_y_sinsum_phase,dof_line)
    reduced_chi_sinsum=cs.reduced_chi_square(frac_rms_sinsum_chi,phase_sinsum_chi,dof_sinsum_rms,dof_sinsum_phase)
    return parameters_sinsum_fracrms,fit_y_sinsum_fracrms,parameters_sinsum_phase,fit_y_sinsum_phase,dof_sinsum_rms,dof_sinsum_phase,frac_rms_sinsum_chi,phase_sinsum_chi,reduced_chi_sinsum

