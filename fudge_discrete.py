import numpy as np 
import load_and_clean as lac
import rms_normalisation as rn
from stingray import Lightcurve, AveragedCrossspectrum
from astropy.io import fits
import filter_fits_gti as ffg
import spin_phase_calculator as phasee
import calc_stokes_analytic as csa
import math
from stingray import Powerspectrum
import matplotlib.pyplot as plt
import os


def fudge_factor(file1,file2,file_1,file_2,file_3,Pmin,Pmax,bin_length,seg_length,fmin,fmax,gti,response_dir,phase_bin_num,v,v_dot,v_double_dot,t_0):
    
    
    #First calculate rms from observed data
    
    
    file_subject_data,file_subject_header,*_=lac.load_and_clean(file1,Pmin,Pmax)
    file_ref_data,file_ref_header,*_=lac.load_and_clean(file2,Pmin,Pmax)

    norm_factor,lc_ref,GTI=rn.rms_normalisation(file_subject_data,file_ref_data,bin_length,seg_length,fmin,fmax,gti)
  
    lc_1_ref=Lightcurve.make_lightcurve(file_subject_data['TIME'],dt=bin_length,gti=GTI)
    lc_1_ref.apply_gtis()
    lc_1_ref_countrate=lc_1_ref.meanrate

    lc_2_ref=Lightcurve.make_lightcurve(file_ref_data['TIME'],dt=bin_length,gti=GTI)
    lc_2_ref.apply_gtis()
    lc_2_ref_countrate=lc_2_ref.meanrate

    cs_ref=AveragedCrossspectrum.from_lightcurve(lc_1_ref,lc_ref,seg_length,norm='frac')

    mask_cs=(cs_ref.freq>fmin) & (cs_ref.freq<fmax)
    power_real_cs=cs_ref.power.real[mask_cs].mean()
    power_im_cs=cs_ref.power.imag[mask_cs].mean()
    #freq_cs=cs_ref.freq[mask_cs]
    #mag_cs=np.sqrt(power_real_cs**2+power_im_cs**2)
    #mean_cs=np.mean(mag_cs)
    frac_rms_observed=np.sqrt(power_real_cs*(fmax-fmin))  #NEW  rms_s*rms_r
    #print('test:',test)
    #frac_rms_observed=mag_cs*norm_factor
    #frac_rms_observed=test
    print('Fractional RMS observed (fundamental): ',frac_rms_observed)


    #Now calculate the rms from the phase folded lightcurve



    #Loacing and cleaning data
    data_1,header_1,*_=lac.load_and_clean(file_1,Pmin,Pmax)
    data_2,header_2,*_=lac.load_and_clean(file_2,Pmin,Pmax)
    data_3,header_3,*_=lac.load_and_clean(file_3,Pmin,Pmax)


    TIMES_1,W_MOM_1,PI_1,q_1,u_1,PHI_1=ffg.filter_fits_gti(data_1,gti)
    TIMES_2,W_MOM_2,PI_2,q_2,u_2,PHI_2=ffg.filter_fits_gti(data_2,gti)
    TIMES_3,W_MOM_3,PI_3,q_3,u_3,PHI_3=ffg.filter_fits_gti(data_3,gti)

    global Aeff_E_DU1, Aeff_E_DU2, Aeff_E_DU3
    global Aeff_mu_E_DU1, Aeff_mu_E_DU2, Aeff_mu_E_DU3
    
    



    #Load in response files

    #DU1

    #MRF
    with fits.open(response_dir+'/mrf/ixpe_d1_20170101_05.mrf') as hdu:
        SPECRESP_MRF_DU1=hdu[1].data
        Aeff_mu_E_DU1=SPECRESP_MRF_DU1['SPECRESP']
    #ARF
    with fits.open(response_dir+'/arf/ixpe_d1_20170101_05.arf') as hdu:
        SPECRESP_ARF_DU1=hdu[1].data
        Aeff_E_DU1=SPECRESP_ARF_DU1['SPECRESP']

    #DU2

    #MRF
    with fits.open(response_dir+'/mrf/ixpe_d2_20170101_05.mrf') as hdu:
        SPECRESP_MRF_DU2=hdu[1].data
        Aeff_mu_E_DU2=SPECRESP_MRF_DU2['SPECRESP']
    #ARF
    with fits.open(response_dir+'/arf/ixpe_d2_20170101_05.arf') as hdu:
        SPECRESP_ARF_DU2=hdu[1].data
        Aeff_E_DU2=SPECRESP_ARF_DU2['SPECRESP']

    #DU3

    #MRF
    with fits.open(response_dir+'/mrf/ixpe_d3_20170101_05.mrf') as hdu:
        SPECRESP_MRF_DU3=hdu[1].data
        Aeff_mu_E_DU3=SPECRESP_MRF_DU3['SPECRESP']

    #ARF
    with fits.open(response_dir+'/arf/ixpe_d3_20170101_05.arf') as hdu:
        SPECRESP_ARF_DU3=hdu[1].data
        Aeff_E_DU3=SPECRESP_ARF_DU3['SPECRESP']

    
    phase_bin_number = phase_bin_num # Number of pulse phase bins
    phase_minimum = 0  # Minimum normalised phase
    phase_maximum = 1  # Maximum normalised phase

    bspace = np.linspace((phase_minimum), (phase_maximum), phase_bin_number + 1)
    phase_bin_list = [(bspace[i - 1], bspace[i]) for i in range(1, len(bspace))]

    

    phase_plot=[np.mean(i) for i in phase_bin_list] #phase bin centers
    av_phase_err=[(i[1]-i[0])/2 for i  in phase_bin_list]

    phase_DU1=phasee.phase_v_double_dot(v,v_dot,v_double_dot,TIMES_1,t_0)
    
    phase_DU2=phasee.phase_v_double_dot(v,v_dot,v_double_dot,TIMES_2,t_0)
    
    phase_DU3=phasee.phase_v_double_dot(v,v_dot,v_double_dot,TIMES_3,t_0)

    
    PD_array=[]
    dPD_array=[]
    
    PA_array=[]
    dPA_array=[]
    
    C_gamma_array=[]
    
    C_gamma_array_12=[]
    C_gamma_array_3=[]


   
    phase_plot=[np.mean(i) for i in phase_bin_list] #phase bin centers
    
    av_phase_err=[(i[1]-i[0])/2 for i  in phase_bin_list]
 
    
   
    # Storage arrays for results
    C_gamma_array, C_gamma_array_12, C_gamma_array_3 = [], [], []
    PD_array, dPD_array, PA_array, dPA_array = [], [], [], []

    
    # Loop over phase bins
    for i in range(phase_bin_num):
        phase_min, phase_max = phase_bin_list[i][0], phase_bin_list[i][1]
   
        # Initialize cumulative variables
        I_NEFF_TOT, dI_NEFF_TOT_sqrd = 0, 0
        Q_NEFF_TOT, dQ_NEFF_TOT_sqrd = 0, 0
        U_NEFF_TOT, dU_NEFF_TOT_sqrd = 0, 0

        for j in range(1, 4):  # Loop over DU1, DU2, DU3
            phase_DU = eval(f'phase_DU{j}')
            phase_cut = (phase_min <= phase_DU) & (phase_DU <= phase_max)

            # Apply mask to select data
            W_MOM_phase_cut = np.ones(np.sum(phase_cut))  # Weight = 1 for all selected events
            PI_phase_cut = eval(f'PI_{j}')[phase_cut]
            q_phase_cut = eval(f'q_{j}')[phase_cut]
            u_phase_cut = eval(f'u_{j}')[phase_cut]
            Aeff_event = np.array([eval(f'Aeff_E_DU{j}')[k - 1] for k in PI_phase_cut])
            Aeff_mu_event = np.array([eval(f'Aeff_mu_E_DU{j}')[k - 1] for k in PI_phase_cut])

            # Compute Stokes parameters
            I_NEFF, dI_NEFF, Q_NEFF, dQ_NEFF, U_NEFF, dU_NEFF = csa.calculate_stokes(
                W_MOM_phase_cut, q_phase_cut, u_phase_cut, Aeff_event, Aeff_mu_event
            )

            # Accumulate for total
            I_NEFF_TOT += I_NEFF
            dI_NEFF_TOT_sqrd += dI_NEFF**2
            Q_NEFF_TOT += Q_NEFF
            dQ_NEFF_TOT_sqrd += dQ_NEFF**2
            U_NEFF_TOT += U_NEFF
            dU_NEFF_TOT_sqrd += dU_NEFF**2

            if j == 1 or j == 2:  # Store sum of DU1 and DU2
                if j == 1:
                    C_gamma_12 = I_NEFF
                else:
                    C_gamma_12 += I_NEFF
            elif j == 3:  # Store DU3 separately
                C_gamma_3 = I_NEFF

        # Compute total errors
        dI_NEFF_TOT = np.sqrt(dI_NEFF_TOT_sqrd)
        dQ_NEFF_TOT = np.sqrt(dQ_NEFF_TOT_sqrd)
        dU_NEFF_TOT = np.sqrt(dU_NEFF_TOT_sqrd)

        # Compute normalized Stokes parameters
        q_phase_cut = Q_NEFF_TOT / I_NEFF_TOT
        u_phase_cut = U_NEFF_TOT / I_NEFF_TOT
        dq_phase_cut = dQ_NEFF_TOT / I_NEFF_TOT
        du_phase_cut = dU_NEFF_TOT / I_NEFF_TOT

        # Compute polarization degree and angle
        PD_phase_cut = np.sqrt(q_phase_cut**2 + u_phase_cut**2)
        dPD_phase_cut = np.sqrt((q_phase_cut * dq_phase_cut)**2 + (u_phase_cut * du_phase_cut)**2) / PD_phase_cut
        PA_phase_cut = 0.5 * math.atan2(u_phase_cut, q_phase_cut)
        dPA_phase_cut = 0.5 * np.sqrt((u_phase_cut * dq_phase_cut)**2 + (q_phase_cut * du_phase_cut)**2) / (PD_phase_cut**2)

        # Store results
        C_gamma_array.append(I_NEFF_TOT)
        C_gamma_array_12.append(C_gamma_12)
        C_gamma_array_3.append(C_gamma_3)
        PD_array.append(PD_phase_cut)
        dPD_array.append(dPD_phase_cut)
        PA_array.append(PA_phase_cut)
        dPA_array.append(dPA_phase_cut)

 
    #T#=1/v #pulse period
    #ulse_phase_time_convert=[pp*T for pp in phase_plot] #Converting pulse phase to time
    

    #gti=[[0,9999999999999999999]]
    phase_fold_lc=Lightcurve(phase_plot,C_gamma_array)#,gti=gti)
    plt.figure()
    plt.title('Folded Lightcurve')
    plt.plot(phase_fold_lc.time,phase_fold_lc.counts,'.')
    plt.show()

    #plt.figure()
    #plt.plot(phase_fold_lc.time,phase_fold_lc.counts,'.')
    #plt.plot(pulse_phase_time_convert,phase_fold_lc.counts,'.')
    #plt.show()
    #from stingray import AveragedPowerspectrum
    ps_phase_fold=Powerspectrum.from_lightcurve(phase_fold_lc,norm='frac')

    plt.figure()
    plt.title('Power Spectrum of phase folded lc')
    plt.plot(ps_phase_fold.freq,ps_phase_fold.power.real,'.')
    plt.xscale('log')
    plt.yscale('log')
    plt.show()


    #ps_fundamental_power=ps_phase_fold.power.real[0]
    #print('ps_fundamental_power',ps_fundamental_power)

    mean_folded_lc=np.mean(phase_fold_lc.counts)
    standard_deviation_folded_lc=np.std(phase_fold_lc.counts)
    frac_rms_folded_lc=standard_deviation_folded_lc/mean_folded_lc
    print('frac rms of folded lc over all harmonics',frac_rms_folded_lc)

    phase_fold_mag=np.sqrt((ps_phase_fold.power.real)**2+(ps_phase_fold.power.imag)**2)


    all_power=np.sum(phase_fold_mag)
    sqrt_all_power=np.sqrt(all_power)
    print('frac rms over all harmonics:',sqrt_all_power)


    fundamental_mag=phase_fold_mag[0]
    sqrt_fundamental_mag=np.sqrt(fundamental_mag)
    harmonic_rms=sqrt_fundamental_mag
    print('fundamental rms from phase folded lc:',sqrt_fundamental_mag)

    fudge_factor=frac_rms_observed/(harmonic_rms)
    print('frac rms observed',frac_rms_observed)
    print('harmonic rms',harmonic_rms)
    #print('harmonic rms',harmonic_rms)
    print('fudge factor:',fudge_factor*1)

    return fudge_factor
