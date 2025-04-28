
import numpy as np
import math
def simulate_ixpe_observations(num_iter,T,qpo_freq,rms,total_count_rate,Q,delta,p,psi,mu,gamma,k,mod_min,mod_max,
                              A_r,B_r,C_r,D_r,E_r,A_p,B_p,C_p,output_file):
    
    r = (1/3) * totcal_count_rate #Referencr band count rate
    
    s = (2/3) * total_count_rate #Subject band count rate
    
    Prn=2/r #Reference band Poisson noise
    
    delta_phi= ()
    
    for i in mod_angle_list:

        mod_min=i[0]
        mod_max=i[1]

        mod_min_rad=np.radians(mod_min)
        mod_max_rad=np.radians(mod_max)
        phi_rad=(mod_min+mod_max)/2
        phi_arr_rad.append(phi_rad)

        phi=(mod_min+mod_max)/2 # Midpoint of modulation angle bin
        phi_arr.append(phi)

        s_phi=(s/k)*(1+(mu*p*np.cos(2*(psi-phi_rad)))) # Subject band count rate (ie total subject band on all mod angles * modulation function)
        s_phi_arr.append(s_phi)

        Psn=2/s_phi # Power spectrum subject band poisson noise (assuming no dead time)
        Psn_arr.append(Psn)

        rms_phi_0=Model_r(phi)
        rms_model.append(rms_phi_0)

        phase_phi=Model_p(phi)
        phase_model.append(phase_phi)
        
    phi_arr=np.array(phi_arr)
    xerr=[4.5]*len(phi_arr)

  


    rms_0_av=np.array(np.mean(rms_model))
    phase_0_av=np.mean(phase_model)


    rms_phi_final=(rms/rms_0_av)*np.array(rms_model)

    Ps=np.array((rms_phi_final**2)/delta) #Subject band Power spectrrum
 
    Pr=np.array((rms**2)/delta)     #Reference band power spectrum with no noise (kindof equivelent to cospec but we assume no dead time here)
  

    mag_G_sqrd=gamma*Pr*Ps

    Pr_t=Pr+Prn #Power spectrum of reference band inc poisson noise



    Ps_t=Ps+Psn_arr #Power spectrum of subject band inc poisson noise


    dG=np.sqrt(((Pr_t)/(2*N))* (f-(mag_G_sqrd/Pr)) )


    phase_err=(dG/np.sqrt(mag_G_sqrd)[0])/(2*np.pi)
    
    
    #Now for the simulation...
    
    for i in range(number_of_iterations):     

        
        ### Simulating rms and phase data points ###
        frac_rms_sim=[]
        phase_sim=[]

        for i in range(len(phi_arr)):
            sim_frac_rms=float(np.random.normal(rms_model[i],dG[i],1))
            frac_rms_sim.append(sim_frac_rms)
            sim_phase=float(np.random.normal(phase_model[i],phase_err[i],1))
            phase_sim.append(sim_phase)

        frac_rms_sim_mega.append(frac_rms_sim)
        phase_sim_mega.append(phase_sim)


        # Fit sinesum model to simulated data

            #Assigning observed data sineusoid parameters as initial guess

        rms_initial_guess_params=np.array([p1_r,p3_r,p4_r,p6_r,p7_r])
        phase_initial_guess_params=np.array([p1_p,p3_p,p4_p,p6_p,p7_p])

        #finding best fit parameters by minizing fit
        fit_rms_params,cov_rms = curve_fit(Model_r_fit, phi_arr, frac_rms_sim, rms_initial_guess_params,maxfev=100000) #fitting free sine to rms
        fit_phase_params,cov_phase = curve_fit(Model_p_fit, phi_arr, phase_sim, phase_initial_guess_params,maxfev=100000) #fitting free sine to rms
        #print('done')

        p1_rf= fit_rms_params[0]
        #p2_rf= fit_rms_params[1] 
        p3_rf =fit_rms_params[1]
        p4_rf= fit_rms_params[2]
        #p5_rf= fit_rms_params[4]
        p6_rf=fit_rms_params[3]
        p7_rf=fit_rms_params[4]

        p1_pf= fit_phase_params[0]
        #p2_pf= fit_phase_params[1] 
        p3_pf =fit_phase_params[1]
        p4_pf= fit_phase_params[2]
        #p5_pf= fit_phase_params[4]
        p6_pf=fit_phase_params[3]
        p7_pf=fit_phase_params[4]


        fit_rms_f = Model_r_fit(phi_arr,p1_rf,p3_rf,p4_rf,p6_rf,p7_rf)
        fit_phase_f=Model_p_fit(phi_arr,p1_pf,p3_pf,p4_pf,p6_pf,p7_pf) #spitting out fit associated datapoints 

        frac_rms_fit_mega.append(fit_rms_f)
        phase_fit_mega.append(fit_phase_f)

        rms_const_initial_guess=[rms]
        phase_const_initial_guess=[phase_0_av]

        const_fit_rms_params,cov_rms = curve_fit(line, phi, frac_rms_sim, rms_const_initial_guess,maxfev=5000) #fitting free sine to rms
        const_fit_phase_params,cov_phase = curve_fit(line, phi, phase_sim, phase_const_initial_guess,maxfev=5000) #fitting free sine to rms

        const_rms_fit=line(phi,const_fit_rms_params)
        const_phase_fit=line(phi,const_fit_phase_params)

        frac_rms_fit_mega_const.append(const_rms_fit)
        phase_fit_mega_const.append(const_phase_fit)



        #stats
        dof_model=len(phi_arr)-5
        dof_line=len(phi_arr)-1

        chi_line_rms = np.sum( ( ( (frac_rms_sim - const_rms_fit)**2 ) / np.array(dG)**2) )
        chi_line_phase = np.sum( ( ( (phase_sim - const_phase_fit)**2 ) / phase_err**2) )

        chi_sine_rms_f=np.sum((frac_rms_sim - fit_rms_f)**2/np.array(dG)**2)
        chi_sine_phase_f=np.sum((phase_sim - fit_phase_f)**2/np.array(dG)**2)




        sig=F_test(chi_line_rms,chi_line_phase,dof_line,chi_sine_rms_f,chi_sine_phase_f,dof_model,dof_model)
    #    
        sig_vals.append(sig)



    sig_av=np.sum([d for d in sig_vals])/len(sig_vals)
    print('av_sig',sig_av)


    
    