#Generates the chi square value of a fit from the real data, the fit and the associated error on the real data points 

def chi_square(data,fit,error):
    chi_square = np.sum( ( ( (data - fit)**2)/error**2)  ) #from classic chi sqr test
    return chi_square

def reduced_chi_square(chi_square_rms,chi_square_phase,dof_rms,dof_phase):
    reduced_chi=(chi_square_rms+chi_square_phase)/(dof_rms+dof_phase)
    return reduced_chi