#Discretised version of QPO finder
import sys
sys.path.append('/home/c2032014/py_files/')
import load_and_clean as lac
import rms_phase_calculator_NEW as rpc
import rms_normalisation as rn
import numpy as np
import Ingram_2019_errors as I_19errs
import fit_models as fm
import fit_rms_phase as frp
import F_test as ft
def qpo_finder(file1,file2,output_file,gti,bin_length,seg_length,fmin,fmax,Pmin,Pmax,mod_bin_number,coherence_corrector=True):

    #load data and clean
    data_1,header_1,*_=lac.load_and_clean(file1,Pmin,Pmax)
    data_2,header_2,*_=lac.load_and_clean(file2,Pmin,Pmax)

    #calculate rms normalisation
    norm_factor,lc_ref=rn.rms_normalisation(data_1,data_2,bin_length,seg_length,fmin,fmax,gti)

    #Calculate G, subtracting spurious polarisation

    av_mod,av_mod_err,frac_rms,frac_rms_err,phase, phase_err=rpc.rms_phase_anderr_calc(data_1,data_2,gti,lc_ref,bin_length,seg_length,fmin,fmax,norm_factor,mod_bin_number,coherence_corrector)
    
    #Statistical significance of fits

    #Fitting null to rms and phase and calculating chi squared values
    parameters_line_fracrms, fit_y_line_fracrms, parameters_line_phase, fit_y_line_phase, dof_line, frac_rms_line_chi, phase_line_chi, reduced_chi_line = frp.fit_line(av_mod, frac_rms, frac_rms_err, phase, phase_err)

    #Similarly for sine 90
    parameters_sin90_fracrms, fit_y_sin90_fracrms, parameters_sin90_phase, fit_y_sin90_phase, dof_sin90_rms, dof_sin90_phase, frac_rms_sin90_chi, phase_sin90_chi, reduced_chi_sin90 = frp.fit_sine_90(av_mod, frac_rms, frac_rms_err, phase, phase_err)

    #Similarly for sine 180
    parameters_sin180_fracrms, fit_y_sin180_fracrms, parameters_sin180_phase, fit_y_sin180_phase, dof_sin180_rms, dof_sin180_phase, frac_rms_sin180_chi, phase_sin180_chi, reduced_chi_sin180 = frp.fit_sine_180(av_mod, frac_rms, frac_rms_err, phase, phase_err)



    #Performing F-test to determine confidence of fits
    #     F test with constant null and Rms:180 Phase: 180
    F_180_all=F_test(chi_line_rms,chi_line_phase,dof_line,chi_sin180_rms,chi_sin180_phase,dof_sin180_rms,dof_sin180_phase)
    print('RMS 180 Phase 180 Significance:',F_180_all)
    print('')

    # F test with constant null and Rms:180 Phase: 90
    F_Rms180_Phase90=F_test(chi_line_rms,chi_line_phase,dof_line,chi_sin180_rms,chi_sin90_phase,dof_sin180_rms,dof_sin90_phase)
    print('RMS 180 Phase 90 Significance:',F_Rms180_Phase90)
    print('')

    # F test with constant null and Rms:180 Phase: sinsum
  
    F_Rms180_Phasesinsum=F_test(chi_line_rms,chi_line_phase,dof_line,chi_sin180_rms,chi_sinsum_phase,dof_sin180_rms,dof_sinsum_phase)
    print('RMS 180 Phase sinsum Significance:',F_Rms180_Phasesinsum)
    print('')


                    # F test with constant null and Rms:90 Phase:180
   
    F_Rms90_Phase180=F_test(chi_line_rms,chi_line_phase,dof_line,chi_sin90_rms,chi_sin180_phase,dof_sin90_rms,dof_sin180_phase)
    print('RMS 90 Phase 180 Significance:',F_Rms90_Phase180)
    print('')

                    # F test with constant null and Rms:90 Phase:90

    F_Rms90_Phase90=F_test(chi_line_rms,chi_line_phase,dof_line,chi_sin90_rms,chi_sin90_phase,dof_sin90_rms,dof_sin90_phase)
    print('RMS 90 Phase 90 Significance:',F_Rms90_Phase90)
    print('')

                    #F test for 90 rms sinsum phase vs constant

    F_90_rms_phase_sinsum=F_test(chi_line_rms,chi_line_phase,dof_line,chi_sin90_rms,chi_sinsum_phase,dof_sin90_rms,dof_sinsum_phase)
    print('RMS 90 Phase sinsum Significance:',F_90_rms_phase_sinsum)
    print('')
                             
                     #F test for sinsum rms phase 180 vs constant


    F_sinsum_rms_phase_180=F_test(chi_line_rms,chi_line_phase,dof_line,chi_sinsum_rms,chi_sin180_phase,dof_sinsum_rms,dof_sin180_phase)
    print('RMS sinsum Phase 180 Significance:',F_sinsum_rms_phase_180)
    print('')

                      
                    #F test for sinsum rms phase 90 vs constant

    F_sinsum_rms_phase_90=F_test(chi_line_rms,chi_line_phase,dof_line,chi_sinsum_rms,chi_sin90_phase,dof_sinsum_rms,dof_sin90_phase)
    print('RMS Sinsum Phase 90 Significance:',F_sinsum_rms_phase_90)
    print('')

                      
                    # Ftest for sinsum rms+phase and constant

    F_sinsum=F_test(chi_line_rms,chi_line_phase,dof_line,chi_sinsum_rms,chi_sinsum_phase,dof_sinsum_rms,dof_sinsum_phase)
    print('RMS sinsum Phase sinsum Significance:',F_sinsum)
    print('')

    results=np.array(tuple(zip(av_mod,av_mod_err,frac_rms,frac_rms_err,phase,phase_err,fit_y_line_fracrms,fit_y_line_phase,fit_y_sin90_rms,fit_y_sin90_phase,fit_y_sin180_rms,fit_y_sin180_phase,fit_y_sinsum_rms,fit_y_sinsum_phase)))
    np.savetxt(output_file,results)
    return results



