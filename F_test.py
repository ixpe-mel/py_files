

#Computes the statistical significance of a model over a given null hypothesis form their chi square values and degrees of freedom
import scipy


def F_test(chi_null_rms,chi_null_phase,dof_null_rms,dof_null_phase,chi_model_rms,chi_model_phase,dof_model_rms,dof_model_phase):

        F= (((chi_null_rms+chi_null_phase)-(chi_model_rms+chi_model_phase))/((dof_null_rms+dof_null_phase)-((dof_model_rms+dof_model_phase)))) / ((chi_model_rms+chi_null_phase)/(dof_null_rms+dof_null_phase))
        print('F',F)
        p_value=1-scipy.stats.f.cdf(F,(dof_null_rms+dof_null_phase),((dof_model_rms+dof_model_phase)))
        print(' P value',p_value)
        confidence  = 1.0 - p_value
        print('Confidence',confidence)
        significance  = 2.0**0.5 * scipy.special.erfinv(confidence)
        print('Significance:',significance)
        return significance
 