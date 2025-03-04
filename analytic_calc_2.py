import numpy as np    
import sys
sys.path.append('/home/c2032014/py_files')
import load_and_clean as lac
import filter_fits_gti as ffg
import importlib
importlib.reload(ffg)
from astropy.io import fits
import spin_phase_calculator as phasee
import calc_stokes_analytic as csa
import matplotlib.pyplot as plt
from stingray import Lightcurve, Powerspectrum, Crossspectrum, AveragedCrossspectrum
import math
importlib.reload(csa)

def Analytic_model_rms_phase(source_name,file_1,file_2,
           file_3,ixpe_master_GTI,
           response_dir,v,v_dot,v_double_dot,t_0,fmin,fmax,Pmin,Pmax,phase_bin_num,mueff):
    
    GTI_master=list(np.loadtxt(str(ixpe_master_GTI)))

    #Loacing and cleaning data
    data_1,header_1,*_=lac.load_and_clean(file_1,Pmin,Pmax)
    data_2,header_2,*_=lac.load_and_clean(file_2,Pmin,Pmax)
    data_3,header_3,*_=lac.load_and_clean(file_3,Pmin,Pmax)


    TIMES_1,W_MOM_1,PI_1,q_1,u_1,PHI_1=ffg.filter_fits_gti(data_1,ixpe_master_GTI)
    TIMES_2,W_MOM_2,PI_2,q_2,u_2,PHI_2=ffg.filter_fits_gti(data_2,ixpe_master_GTI)
    TIMES_3,W_MOM_3,PI_3,q_3,u_3,PHI_3=ffg.filter_fits_gti(data_3,ixpe_master_GTI)

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
        I_TOT, dI_TOT_sqrd = 0, 0
        Q_TOT, dQ_TOT_sqrd = 0, 0
        U_TOT, dU_TOT_sqrd = 0, 0

        Q_uncorr, U_uncorr, I_uncorr = 0, 0, 0

        for j in range(1, 4):  # Loop over DU1, DU2, DU3
            phase_DU = eval(f'phase_DU{j}')
            phase_cut = (phase_min <= phase_DU) & (phase_DU <= phase_max)

            # Apply mask to select data
            W_MOM_phase_cut = np.ones(np.sum(phase_cut))  # Weight = 1 for all selected events#
            #W_MOM_phase_cut=eval(f'W_MOM_{j}')[phase_cut]
            #print(len(W_MOM_phase_cut))
            PI_phase_cut = eval(f'PI_{j}')[phase_cut]
            q_phase_cut = eval(f'q_{j}')[phase_cut]
            u_phase_cut = eval(f'u_{j}')[phase_cut]
            Aeff_event = np.array([eval(f'Aeff_E_DU{j}')[k - 1] for k in PI_phase_cut])
            Aeff_mu_event = np.array([eval(f'Aeff_mu_E_DU{j}')[k - 1] for k in PI_phase_cut])

            #Uncorrected (for modulation fa ctor) stokes

            q_phase_cut_uncorrected=np.sum(q_phase_cut)
            u_phase_cut_uncorrected=np.sum(u_phase_cut)
            i_phase_cut_uncorrected=np.sum(W_MOM_phase_cut)

            Q_uncorr+=q_phase_cut_uncorrected
            U_uncorr+=u_phase_cut_uncorrected
            I_uncorr+=i_phase_cut_uncorrected

            # Compute Stokes parameters
            I, dI, Q, dQ, U, dU = csa.calculate_stokes(
                W_MOM_phase_cut, q_phase_cut, u_phase_cut, Aeff_event, Aeff_mu_event
            )

            # Accumulate for total
            I_TOT += I
            dI_TOT_sqrd += dI**2
            Q_TOT += Q
            dQ_TOT_sqrd += dQ**2
            U_TOT += U
            dU_TOT_sqrd += dU**2

            if j == 1 or j == 2:  # Store sum of DU1 and DU2
                if j == 1:
                    C_gamma_12 = I
                else:
                    C_gamma_12 += I
            elif j == 3:  # Store DU3 separately
                C_gamma_3 = I

        # Compute total errors
        dI_TOT = np.sqrt(dI_TOT_sqrd)
        dQ_TOT = np.sqrt(dQ_TOT_sqrd)
        dU_TOT = np.sqrt(dU_TOT_sqrd)

        # Compute normalized Stokes parameters
        q_phase_cut = Q_TOT / I_TOT
        u_phase_cut = U_TOT / I_TOT
        dq_phase_cut = dQ_TOT / I_TOT
        du_phase_cut = dU_TOT / I_TOT

        # Compute polarization degree and angle
        PD_phase_cut = np.sqrt(q_phase_cut**2 + u_phase_cut**2)
        dPD_phase_cut = np.sqrt((q_phase_cut * dq_phase_cut)**2 + (u_phase_cut * du_phase_cut)**2) / PD_phase_cut
        PA_phase_cut = 0.5 * math.atan2(u_phase_cut, q_phase_cut)
        dPA_phase_cut = 0.5 * np.sqrt((u_phase_cut * dq_phase_cut)**2 + (q_phase_cut * du_phase_cut)**2) / (PD_phase_cut**2)

        # Store results
        C_gamma_array.append(I_TOT)
        C_gamma_array_12.append(C_gamma_12)
        C_gamma_array_3.append(C_gamma_3)
        PD_array.append(PD_phase_cut)
        dPD_array.append(dPD_phase_cut)
        PA_array.append(PA_phase_cut)
        dPA_array.append(dPA_phase_cut)

    PA_array=np.array(PA_array)
    PD_array=np.array(PD_array)


    #Calculating effective modulation factor


    #mu_eff_calculated=np.sqrt((Q_uncorr/I_uncorr)**2 + (U_uncorr/I_uncorr)**2)/np.mean(PD_array)
    #mu_eff=np.sqrt(  (Q_uncorr**2+U_uncorr**2)/ (Q_TOT**2+U_TOT**2)    )
    #print('mu_eff:',mu_eff_calculated)

    
    











    #T=1/v #pulse period
    #pulse_phase_time_convert=[pp*T for pp in phase_plot] #Converting pulse phase to time
    
    #plt.figure()
    #plt.title('I')
    #plt.plot(pulse_phase_time_convert,C_gamma_array,'.')
    #plt.show()

    #plt.figure()
    #plt.plot(phase_plot,C_gamma_array,'.')
    #plt.show()

    #plt.figure()
    #plt.title('PD')
    #plt.plot(pulse_phase_time_convert,PD_array,'.')
    #plt.show()

    #plt.figure()
    #plt.title('PA')
    #plt.plot(pulse_phase_time_convert,PA_array,'.')
    #plt.show()
    #
    C_gamma_array_orig=C_gamma_array
    lc_I_orig=Lightcurve(phase_plot,C_gamma_array_orig)
    I_ps=Powerspectrum.from_lightcurve(lc_I_orig,norm='frac')
    I_ps_freqcut=I_ps.power.real[0]
    I_ps_frac_rms=np.sqrt(I_ps_freqcut)
    print('Frac rms of C(\gamma) fundamental:',I_ps_frac_rms)


     #New fudge 
    rms_ratio=0.28520197995386165/I_ps_frac_rms #rxj estimate
    print('rxj rms ratio:',rms_ratio)
    #rms_ratio=1.2847  #herx1 estimate
    print('rms ratio:',rms_ratio)
    
    C_gamma_array=((C_gamma_array-np.mean(C_gamma_array))*rms_ratio)+np.mean(C_gamma_array)

    plt.figure()
    plt.title('I scaling')
    plt.plot(phase_plot,C_gamma_array,'.',label='new')
    plt.plot(phase_plot,C_gamma_array_orig,'.',label='old')
    plt.legend()
    plt.show()


    #lc_ref_3=Lightcurve(phase_plot,C_gamma_array_3)          
    #lc_ref_12=Lightcurve(phase_plot,C_gamma_array_12)       
   
    #lc_ref_3_mean=np.mean(C_gamma_array_3)
    #lc_ref_3_std=np.std(C_gamma_array_3)
    #frac_rms_ref_3=lc_ref_3_std/lc_ref_3_mean
    #print('frac rms of DU3 reference lc:',frac_rms_ref_3)



    lc_I=Lightcurve(phase_plot,C_gamma_array)

    #lc_ref_3_ps=Powerspectrum.from_lightcurve(lc_ref_3,norm='frac')
    #lc_ref_3_ps=Powerspectrum.from_lightcurve(lc_ref_3,norm='frac')
    #print('lc_ref_3_ps',lc_ref_3_ps.power.real[:10])
    #lc_ref_3_ps_freqcut=lc_ref_3_ps.power.real[0]

    #print('lc_ref_3_ps 0',lc_ref_3_ps_freqcut)

    #cs_test=Crossspectrum.from_lightcurve(lc_ref_3,lc_ref_3,norm='frac')
    #print('cs_0',cs_test.power.real[0])

    ps_ref=Powerspectrum.from_lightcurve(lc_I,norm='frac')
    ps_ref_freqcut=ps_ref.power.real[0]
    ps_ref_imag=ps_ref.power.imag[0]
    ps_ref_mag=np.sqrt(ps_ref_freqcut**2 + ps_ref_imag**2)

    norm=1/np.sqrt(ps_ref_mag)


    #cs_ref=Crossspectrum.from_lightcurve(lc_ref_12,lc_ref_3,norm='frac')
    #cs_ref_real=cs_ref.power.real[0]
    #cs_ref_imag=cs_ref.power.imag[0]
    #cs_ref_mag=np.sqrt(cs_ref_real**2 + cs_ref_imag**2)
    
    #SWITCH 2 PS
    
    #print('cs_ref_real',cs_ref_real[:10])
    #cs_ref_norm=cs_ref_real[0]
    #print('cs',cs_ref_norm)
    #plt.figure()
    #plt.plot(cs_ref.freq,cs_ref_real,'.')
    #plt.title('cs_ref_real')
    #plt.xscale('log')
    #plt.yscale('log')
    #plt.show()

  
 #   norm=1/(np.sqrt(cs_ref_mag))
 #   print('norm',norm)

   


    #Sanity check 1: frac rms calculted from I
    mean_I=np.mean(C_gamma_array)
    std_I=np.std(C_gamma_array)
    frac_rms_I=std_I/mean_I 
    print('Frac rms of C(\gamma):',frac_rms_I)


    mean_I_12=np.mean(C_gamma_array_12)
    std_I_12=np.std(C_gamma_array_12)
    frac_rms_I_12=std_I_12/mean_I_12
    print('Frac rms of C(\gamma) DU1+DU2:',frac_rms_I_12)

    
    #Sanity check 2: frac rms calculated from ps(I)
    I_ps=Powerspectrum.from_lightcurve(lc_I,norm='frac')
    I_ps_freqcut=I_ps.power.real[0]
    I_ps_frac_rms=np.sqrt(I_ps_freqcut)
    print('Frac rms of C(\gamma) fundamental:',I_ps_frac_rms)

   


    












    #plt.figure()
    #plt.title('Power Spectrum of I phase')
    #plt.plot(I_ps_phase.freq,I_ps_phase.power.real,'.',label='phase')
    #plt.plot(I_ps.freq,I_ps.power.real,'.',label='timec')
    #plt.xscale('log')
    #plt.yscale('log')
    #plt.show()


    #Sanity check 2.5: cs

   # I_cs=Crossspectrum.from_lightcurve(lc_ref_12,lc_ref_3,norm='frac')
   # I_cs_mag=[np.sqrt(i.real**2 + i.imag**2) for i in I_cs.power]
   # frac_rms_cs_I=np.sqrt(np.sum(I_cs_mag))
   # print('I cs all harmonics frac rms:',frac_rms_cs_I)
   # I_cs_freqcut=I_cs_mag[0]
    #[(fmin<=I_cs.freq) & (I_cs.freq<=fmax)]

   # I_cs_frac_rms=np.sqrt(I_cs_freqcut)
   # print('I cs fundamental frac rms',I_cs_frac_rms)



    









    mod_angle_bin_number = 20  # Number of modulation angle bins
    mod_angle_minimum = np.radians(-90)  # Minimum modulation angle in degrees
    mod_angle_maximum = np.radians(90)  # Maximum modulation angle in degrees

    modspace = np.linspace((mod_angle_minimum), (mod_angle_maximum), mod_angle_bin_number + 1)
    mod_angle_bin_list = [(modspace[i - 1], modspace[i]) for i in range(1, len(modspace))]
    mod_angle_bin_width=mod_angle_bin_list[1][1]-mod_angle_bin_list[1][0]
    mod_angle_bin_midpoint=np.array([np.mean(i) for i in mod_angle_bin_list])
    

    mod_angle_bin_list=enumerate(mod_angle_bin_list,1)
    #print('mod_angle_bin_list',mod_angle_bin_list)

    cs_12_3_real_arr=[]
    cs_12_3_imag_arr=[]
    cs_12_3_mod=[]
    frac_rms_dc_12_arr=[]
    frac_rms_dc_12_ps_arr=[]
    mod_angle_bin_midpoint=[]

    #plt.figure()
    #plt.plot(phase_plot,C_gamma_array_3,'.')
    #plt.plot(lc_ref_3.time,lc_ref_3.counts,'.')
    #plt.show()


    plt.figure()
    plt.plot(phase_plot,PD_array,'.')
    plt.title('PD')
    plt.show()

    plt.figure()
    plt.plot(phase_plot,PA_array,'.')
    plt.title('PA')
    plt.show()

    for i in mod_angle_bin_list:

        mod_mid=np.mean(i[1])
        mod_angle_bin_midpoint.append(mod_mid)

        #mod_func= (mod_angle_bin_width/np.pi) * ( 1 + ( (0.3*PD_array) * (np.cos((2*PA_array)-(2*mod_mid))) ))
        #plt.figure()
        #plt.title('mod func')
        #plt.plot(phase_plot,mod_func,'.')
        #plt.show()

        dc= np.array(C_gamma_array) * (mod_angle_bin_width/np.pi) * ( 1 + ( (mueff*PD_array) * (np.cos((2*PA_array)-(2*mod_mid))) ))
        
        
        I_lc_12_subject=Lightcurve(phase_plot,dc)
        #I_lc_12_subject_ps=Powerspectrum.from_lightcurve(I_lc_12_subject,norm='frac')
        #I_lc_12_subject_cs=Crossspectrum.from_lightcurve(I_lc_12_subject,I_lc_12_subject,norm='frac')
        I_lc_12_sub_mean=np.mean(I_lc_12_subject.counts)
        I_lc_12_sub_std=np.std(I_lc_12_subject.counts)
        frac_rms_sub=I_lc_12_sub_std/I_lc_12_sub_mean
        #print('frac rms sub:',frac_rms_sub)

        
        #print('I_lc_12_subject_cs',np.sqrt(I_lc_12_subject_cs.power.real[0]))
        #frac_rms_ps_I_lc_12_subject=np.sqrt(np.sum(I_lc_12_subject_ps.power))
        #print('I_lc_12_subject_ps frac rms all harm:',frac_rms_ps_I_lc_12_subject)
        #frac_rms_dc_12_ps_arr.append(frac_rms_ps_I_lc_12_subject)
        #dc_12_mean=np.mean(dc_12)
        #dc_12_std=np.std(dc_12)
        #frac_rms_dc_12=dc_12_std/dc_12_mean
        #print('dc_12 frac rms:',frac_rms_dc_12)
        #frac_rms_dc_12_arr.append(frac_rms_dc_12)
        
        #plt.figure()
        #plt.title('Analytic subject lc')
        #plt.plot(I_lc_12_subject.time,I_lc_12_subject.counts)
        #plt.show()
    
        cs_12_3=Crossspectrum.from_lightcurve(I_lc_12_subject,lc_I,norm='frac')
        cs_12_3_real=cs_12_3.power.real
        #print('cs_12_3_real',cs_12_3_real[:10])
        cs_12_3_imag=cs_12_3.power.imag
        
        mag_cs=(np.sqrt(i**2+j**2) for i,j in zip(cs_12_3_real,cs_12_3_imag))

        frac_rms_mag_cs=np.sqrt(np.sum(mag_cs))
        frac_rms_cs=np.sqrt(np.sum(cs_12_3.power.real))
        #print('frac rms cs:',frac_rms_cs)   
        #print('frac rms mag cs:',frac_rms_mag_cs)
        #plt.figure()
        #plt.title('cs_12_3_real')
        #plt.plot(cs_12_3.freq,cs_12_3_real,'.')
        #plt.xscale('log')
        #plt.yscale('log')
        #plt.show()

        #plt.figure()
        #plt.title('cs_12_3_imag')
        #plt.plot(cs_12_3.freq,cs_12_3_imag,'.')
        #plt.xscale('log')
        #plt.yscale('log')
        #plt.show()

        cs_12_real_freqcut=cs_12_3_real[0]
        cs_12_imag_freqcut=cs_12_3_imag[0]



       # cs_12_freqcut_mod=np.sqrt(cs_12_real_freqcut**2 + cs_12_imag_freqcut**2)
        #cs_12_3_mod.append(cs_12_freqcut_mod)


        cs_12_3_real_arr.append(cs_12_real_freqcut)
        cs_12_3_imag_arr.append(cs_12_imag_freqcut)

        #cs_12_3_mag=[np.sqrt(i.real**2 + i.imag**2) for i in cs_12_3.power]
        #cs_12_3_freqcut=cs_12_3_mag[0]
        #[(fmin<=cs_12_3_time_convert.freq) & (cs_12_3_time_convert.freq<=fmax)].mean()
        #cs_12_3_arr.append(cs_12_3_freqcut)
        #cs_12_3_freqcut_imag=cs_12_3_time_convert.power.imag[0]
        #[(fmin<=cs_12_3_time_convert.freq) & (cs_12_3_time_convert.freq<=fmax)].mean()
        #cs_12_3_time_convert_imag.append(cs_12_3_freqcut_imag)
       
        #cs_12_3_0.append(complex(cs_12_3.power.real[0],cs_12_3.power.imag[0]))
       
    #print('cs_12_3_time_convert_real',cs_12_3_time_convert_real)
    #cs_12_33=[complex(r, i) for r, i in zip(cs_12_3_time_convert_real, cs_12_3_time_convert_imag)]
    #mod_G_12_3_time_convert=[np.sqrt(i.real**2 + i.imag**2) for i in cs_12_33]
    
    
    

    mod_G=[np.sqrt(i**2 + j**2) for i,j in zip(cs_12_3_real_arr,cs_12_3_imag_arr)]
    print('mod_G',mod_G)
    rms_12_3=[norm*i for i in mod_G]
    

    atan2_vec=np.vectorize(math.atan2)
    phase = [(atan2_vec(i,j))/(2*np.pi) for i,j in zip(cs_12_3_imag_arr,cs_12_3_real_arr)]

    #plt.figure()
    #plt.title('frac rms dc12')
    #plt.plot(mod_angle_bin_midpoint,frac_rms_dc_12_arr,'.',label='frac rms dc 12')
    #plt.show()
    #plt.figure()
    #plt.plot(mod_angle_bin_midpoint,rms_12_3,'.',label='rms123')
    #plt.show()

    #plt.figure()
    #plt.title('frac rms dc12 ps')
    #plt.plot(mod_angle_bin_midpoint,frac_rms_dc_12_ps_arr,'.',label='frac rms dc 12 ps')
    #plt.show()
    #plt.show()

    plt.figure()
    plt.title('rms')
    plt.plot(mod_angle_bin_midpoint,rms_12_3,'.')
    plt.show()
    plt.figure()
    plt.title('phase')
    plt.plot(mod_angle_bin_midpoint,phase,'.')
    plt.show()
     
    #print('rms',rms)
    #print('phase',phase)
        
        
    #Converting phase bin PA/PD to degrees and percent  
        
    PA_array=[np.degrees(i) for i in PA_array]
    dPA_array=[np.degrees(i) for i in dPA_array]
        
    PD_array=[100*i for i in PD_array]
    dPD_array=[100*i for i in dPD_array]
        
    
    
    #analytic_saveout=np.array([rms,phase])
    analytic_saveout=np.array(tuple(zip(phase_plot, av_phase_err, PD_array, dPD_array, PA_array, dPA_array,mod_angle_bin_midpoint,phase,rms_12_3)))
    np.savetxt(f'/home/c2032014/Analytic_calculation/{source_name}_{phase_bin_num}_analytic_calc.txt',analytic_saveout)

    
    return phase_plot, av_phase_err, PD_array, dPD_array, PA_array, dPA_array, mod_angle_bin_midpoint,phase,rms_12_3






        

        