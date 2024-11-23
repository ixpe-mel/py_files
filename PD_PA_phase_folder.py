#PD/PA Folding over event by event phase
import numpy as np
from astropy.io import fits
import math
import matplotlib.pyplot as plt
from stingray import Lightcurve
import importlib

def PD_PA_PHASE(ixpe_event_file_DU1,
           ixpe_event_file_DU2,
           ixpe_event_file_DU3,ixpe_master_GTI,
           response_dir,Pmin,Pmax,phase_bin_number,v,v_dot,v_double_dot,t_0,output_file):
    

    GTI_master=list(np.loadtxt(str(ixpe_master_GTI)))
    
    import sys
    sys.path.append('/home/c2032014/PhD/py_files')
    import spin_phase_calculator as phase
    importlib.reload(phase)
        # Event file DU1
        
    with fits.open(ixpe_event_file_DU1) as hdu:
        EVENTS_DU1=hdu[1].data
        EVENTS_DU1=EVENTS_DU1[(Pmin<=EVENTS_DU1['PI']) & (EVENTS_DU1['PI']<= Pmax)]
        header=hdu[1].header
        TSTART=header['TSTART']
        W_MOM_DU1=list(EVENTS_DU1['W_MOM'])
        PI_DU1=list(EVENTS_DU1['PI'])
        q_1=list(EVENTS_DU1['Q'])
        u_1=list(EVENTS_DU1['U'])
        TIMES_DU1=list(EVENTS_DU1['TIME'])
        
     #Using lc to filter events over gti

        lc_DU1 = Lightcurve(EVENTS_DU1['TIME'], [1]*len(EVENTS_DU1['TIME']) ,gti=GTI_master)
        lc_DU1.apply_gtis()
        filter_set_DU1=set(lc_DU1.time)
    
        events_1=[(index,value) for index,value in enumerate(EVENTS_DU1['TIME']) if value in filter_set_DU1]
        #print(events_1[:100])
        indicies_DU1=list([index for index, value in events_1])
        TIMES_DU1=[value for index,value in events_1]
        print('number of filtered events from DU1',len(TIMES_DU1))
        W_MOM_DU1=np.array([W_MOM_DU1[i] for i in indicies_DU1])
        PI_DU1=np.array([PI_DU1[i] for i in indicies_DU1])
        q_1=np.array([q_1[i] for i in indicies_DU1])
        u_1=np.array([u_1[i] for i in indicies_DU1])


        #events=[x for x in EVENTS_DU1['TIME'] if x in filter_set]
  
    
        # Event file DU2
    
    with fits.open(ixpe_event_file_DU2) as hdu:
        EVENTS_DU2=hdu[1].data
        EVENTS_DU2=EVENTS_DU2[(Pmin<=EVENTS_DU2['PI']) & (EVENTS_DU2['PI']<= Pmax)]
        GTI_DU2=hdu[2].data
        
        W_MOM_DU2=EVENTS_DU2['W_MOM']
        PI_DU2=EVENTS_DU2['PI']
        q_2=EVENTS_DU2['Q']
        u_2=EVENTS_DU2['U']
        
        
        lc_DU2 =Lightcurve(EVENTS_DU2['TIME'], [1]*len(EVENTS_DU2['TIME']) ,gti=GTI_master)
        lc_DU2.apply_gtis()
        
        
      
        filter_set_DU2=set(lc_DU2.time)
    
    
        
        events_2=[(index,value) for index,value in enumerate(EVENTS_DU2['TIME']) if value in filter_set_DU2]
        indicies_DU2=list([index for index, value in events_2])
        TIMES_DU2=[value for index,value in events_2]
        print('number of filtered events from DU2',len(TIMES_DU2))
        W_MOM_DU2=np.array([W_MOM_DU2[i] for i in indicies_DU2])
        PI_DU2=np.array([PI_DU2[i] for i in indicies_DU2])
        q_2=np.array([q_2[i] for i in indicies_DU2])
        u_2=np.array([u_2[i] for i in indicies_DU2])
        

        
        # Event file DU3
        
    with fits.open(ixpe_event_file_DU3) as hdu:
        EVENTS_DU3=hdu[1].data
        EVENTS_DU3=EVENTS_DU3[(Pmin<=EVENTS_DU3['PI']) & (EVENTS_DU3['PI']<= Pmax)]
        GTI_DU3=hdu[2].data
        
        W_MOM_DU3=EVENTS_DU3['W_MOM']
        PI_DU3=EVENTS_DU3['PI']
        q_3=EVENTS_DU3['Q']
        u_3=EVENTS_DU3['U']

        
                
        lc_DU3 =Lightcurve(EVENTS_DU3['TIME'], [1]*len(EVENTS_DU3['TIME']) ,gti=GTI_master)
        lc_DU3.apply_gtis()
        
        
      
        filter_set_DU3=set(lc_DU3.time)
    
    
        
        events_3=[(index,value) for index,value in enumerate(EVENTS_DU3['TIME']) if value in filter_set_DU3]
        indicies_DU3=list([index for index, value in events_3])
        TIMES_DU3=[value for index,value in events_3]
        print('number of filtered events from DU3',len(TIMES_DU3))
        
        W_MOM_DU3=np.array([W_MOM_DU3[i] for i in indicies_DU3])
        PI_DU3=np.array([PI_DU3[i] for i in indicies_DU3])
        q_3=np.array([q_3[i] for i in indicies_DU3])
        u_3=np.array([u_3[i] for i in indicies_DU3])
        
            
  
        
        
        
        #DU1 RESPONSE FILES
    
    
    
    #MRF
    with fits.open(response_dir+'/mrf/ixpe_d1_20170101_alpha075_05.mrf') as hdu:
        SPECRESP_MRF_DU1=hdu[1].data
        Aeff_mu_E_DU1=SPECRESP_MRF_DU1['SPECRESP']
        #print(Aeff_mu_E_DU1[0][0])
        
    #ARF
    with fits.open(response_dir+'/arf/ixpe_d1_20170101_alpha075_05.arf') as hdu:
        SPECRESP_ARF_DU1=hdu[1].data
        Aeff_E_DU1=SPECRESP_ARF_DU1['SPECRESP']
        
        ENERGY_LO=SPECRESP_MRF_DU1['ENERG_LO']
        ENERGY_HI=SPECRESP_MRF_DU1['ENERG_HI']
        
       
    
        #DU2 RESPONSE FILES
    
    #MRF
    with fits.open(response_dir+'/mrf/ixpe_d2_20170101_alpha075_05.mrf') as hdu:
        SPECRESP_MRF_DU2=hdu[1].data
        Aeff_mu_E_DU2=SPECRESP_MRF_DU2['SPECRESP']
        
        
    #ARF
    with fits.open(response_dir+'/arf/ixpe_d2_20170101_alpha075_05.arf') as hdu:
        SPECRESP_ARF_DU2=hdu[1].data
        Aeff_E_DU2=SPECRESP_ARF_DU2['SPECRESP']
        
        
        
            #DU3 RESPONSE FILES
    #MRF
    with fits.open(response_dir+'/mrf/ixpe_d3_20170101_alpha075_05.mrf') as hdu:
        SPECRESP_MRF_DU3=hdu[1].data
        Aeff_mu_E_DU3=SPECRESP_MRF_DU3['SPECRESP']
    #ARF
    with fits.open(response_dir+'/arf/ixpe_d3_20170101_alpha075_05.arf') as hdu:
        SPECRESP_ARF_DU3=hdu[1].data
        Aeff_E_DU3=SPECRESP_ARF_DU3['SPECRESP']
        
        

    #ENERGY_BINS=np.array( list( enumerate( tuple( zip(ENERGY_LO,ENERGY_HI) ),1 ))) 
       
        
    
    
    # Define the parameters

    phase_minimum = 0  # Minimum modulation angle in degrees
    phase_maximum = 2*np.pi  # Maximum modulation angle in degrees

    bspace = np.linspace((phase_minimum), (phase_maximum), phase_bin_number + 1)
    phase_bin_list = [(bspace[i - 1], bspace[i]) for i in range(1, len(bspace))]

    
    phase_plot=[np.mean(i) for i in phase_bin_list]
    
    av_phase_err=[(i[1]-i[0])/2 for i  in phase_bin_list]
    
    #calculating phase for each event for each DU

    
    #RXJ0440 obs 2 (long)
    #v=4.8670*10**-3
    #v_dot=1.42*10**-11
    #v_double_dot=-8*10**-18
    #t_0=1.9393835*10**8
    
    
    
    #Her X-1
    #v=0.8079478799630468
    #v=1/1.2377093
    #v_dot=0
    #v_double_dot=0
    #t_0=TSTART
    
    phase_DU1=phase.phase_v_double_dot(v,v_dot,v_double_dot,np.array(TIMES_DU1),t_0)
    
    phase_DU2=phase.phase_v_double_dot(v,v_dot,v_double_dot,np.array(TIMES_DU2),t_0)
    
    phase_DU3=phase.phase_v_double_dot(v,v_dot,v_double_dot,np.array(TIMES_DU3),t_0)
    
    
    PD_array=[]
    dPD_array=[]
    
    PA_array=[]
    dPA_array=[]

    I_array=[]
    
    #phase_plot=[]
    
    for i in phase_bin_list:
        phase_min=i[0]
        phase_max=i[1]
        
        
       
        # Selecting photons within modulation angle bin range
        W_MOM_DU1_phase_cut=W_MOM_DU1[(phase_min<=phase_DU1) & (phase_DU1<=phase_max)]
        PI_DU1_phase_cut=PI_DU1[(phase_min<=phase_DU1) & (phase_DU1<=phase_max)]
        q_1_phase_cut=q_1[(phase_min<=phase_DU1) & (phase_DU1<=phase_max)]
        u_1_phase_cut=u_1[(phase_min<=phase_DU1) & (phase_DU1<=phase_max)]
                          
                        
        W_MOM_DU2_phase_cut=W_MOM_DU2[(phase_min<=phase_DU2) & (phase_DU2<=phase_max)]
        PI_DU2_phase_cut=PI_DU2[(phase_min<=phase_DU2) & (phase_DU2<=phase_max)]
        q_2_phase_cut=q_2[(phase_min<=phase_DU2) & (phase_DU2<=phase_max)]
        u_2_phase_cut=u_2[(phase_min<=phase_DU2) & (phase_DU2<=phase_max)]
                          
                          
        W_MOM_DU3_phase_cut=W_MOM_DU3[(phase_min<=phase_DU3) & (phase_DU3<=phase_max)]
        PI_DU3_phase_cut=PI_DU3[(phase_min<=phase_DU3) & (phase_DU3<=phase_max)]
        q_3_phase_cut=q_3[(phase_min<=phase_DU3) & (phase_DU3<=phase_max)]
        u_3_phase_cut=u_3[(phase_min<=phase_DU3) & (phase_DU3<=phase_max)]
        
        # Assigning response variables to each event   
        Aeff_mu_event_DU1=[Aeff_mu_E_DU1[i-1] for i in PI_DU1_phase_cut]   
        Aeff_mu_event_DU2=[Aeff_mu_E_DU2[i-1] for i in PI_DU2_phase_cut]
        Aeff_mu_event_DU3=[Aeff_mu_E_DU3[i-1] for i in PI_DU3_phase_cut]

        Aeff_event_DU1=[Aeff_E_DU1[i-1] for i in PI_DU1_phase_cut]
        Aeff_event_DU2=[Aeff_E_DU2[i-1] for i in PI_DU2_phase_cut]
        Aeff_event_DU3=[Aeff_E_DU3[i-1] for i in PI_DU3_phase_cut]
        
        W_MOM_DU1_phase_cut=np.array(list([1]*len(Aeff_event_DU1)))
        W_MOM_DU2_phase_cut=np.array(list([1]*len(Aeff_event_DU2)))
        W_MOM_DU3_phase_cut=np.array(list([1]*len(Aeff_event_DU3)))
        
        #Calculating weighted Stokes
        #DU1
        I_NEFF_1_phase_cut=np.sum(W_MOM_DU1_phase_cut/Aeff_event_DU1)
        dI_NEFF_1_sqrd_phase_cut=np.sum( (W_MOM_DU1_phase_cut/Aeff_event_DU1)**2 ) 
            #DU2
        I_NEFF_2_phase_cut=np.sum(W_MOM_DU2_phase_cut/Aeff_event_DU2)
        dI_NEFF_2_sqrd_phase_cut=np.sum( (W_MOM_DU2_phase_cut/Aeff_event_DU2)**2 ) 
            #DU3
        I_NEFF_3_phase_cut=np.sum(W_MOM_DU3_phase_cut/Aeff_event_DU3)
        dI_NEFF_3_sqrd_phase_cut=np.sum( (W_MOM_DU3_phase_cut/Aeff_event_DU3)**2 ) 
            #TOT
        dI_NEFF_TOT_sqrd_phase_cut=dI_NEFF_1_sqrd_phase_cut+dI_NEFF_2_sqrd_phase_cut+dI_NEFF_3_sqrd_phase_cut
        dI_NEFF_TOT_phase_cut=np.sqrt(dI_NEFF_TOT_sqrd_phase_cut)
        I_NEFF_TOT_phase_cut=I_NEFF_1_phase_cut+I_NEFF_2_phase_cut+I_NEFF_3_phase_cut
        I_array.append(I_NEFF_TOT_phase_cut)



        # Stokes Q

            #DU1
        Q_NEFF_1_phase_cut=np.sum( (W_MOM_DU1_phase_cut*q_1_phase_cut) / Aeff_mu_event_DU1)
        dQ_NEFF_1_sqrd_phase_cut=np.sum( ( (W_MOM_DU1_phase_cut*q_1_phase_cut)/Aeff_event_DU1 )**2 )
            #DU2
        Q_NEFF_2_phase_cut=np.sum( (W_MOM_DU2_phase_cut*q_2_phase_cut) / Aeff_mu_event_DU2)
        dQ_NEFF_2_sqrd_phase_cut=np.sum( ( (W_MOM_DU2_phase_cut*q_2_phase_cut)/Aeff_event_DU2 )**2 )
            #DU3
        Q_NEFF_3_phase_cut=np.sum( (W_MOM_DU3_phase_cut*q_3_phase_cut) / Aeff_mu_event_DU3)
        dQ_NEFF_3_sqrd_phase_cut=np.sum( ( (W_MOM_DU3_phase_cut*q_3_phase_cut)/Aeff_event_DU3 )**2 )
            #TOT
            
            
            
        dQ_NEFF_TOT_sqrd_phase_cut=dQ_NEFF_1_sqrd_phase_cut+dQ_NEFF_2_sqrd_phase_cut+dQ_NEFF_3_sqrd_phase_cut
        dQ_NEFF_TOT_phase_cut=np.sqrt(dQ_NEFF_TOT_sqrd_phase_cut)
        
        
        
        
        
        Q_NEFF_TOT_phase_cut=Q_NEFF_1_phase_cut+Q_NEFF_2_phase_cut+Q_NEFF_3_phase_cut



       # Stokes U 
            #DU1
        U_NEFF_1_phase_cut=np.sum( (W_MOM_DU1_phase_cut*u_1_phase_cut) / Aeff_mu_event_DU1)
        dU_NEFF_1_sqrd_phase_cut=np.sum( ( (W_MOM_DU1_phase_cut*u_1_phase_cut)/Aeff_event_DU1 )**2 )
            #DU2
        U_NEFF_2_phase_cut=np.sum( (W_MOM_DU2_phase_cut*u_2_phase_cut) / Aeff_mu_event_DU2)
        dU_NEFF_2_sqrd_phase_cut=np.sum( ( (W_MOM_DU2_phase_cut*u_2_phase_cut)/Aeff_event_DU2 )**2 )
            #DU3
        U_NEFF_3_phase_cut=np.sum( (W_MOM_DU3_phase_cut*u_3_phase_cut) / Aeff_mu_event_DU3)
        dU_NEFF_3_sqrd_phase_cut=np.sum( ( (W_MOM_DU3_phase_cut*u_3_phase_cut)/Aeff_event_DU3 )**2 )
        
        
        
        
        
        #TOT
        dU_NEFF_TOT_sqrd_phase_cut=dU_NEFF_1_sqrd_phase_cut+dU_NEFF_2_sqrd_phase_cut+dU_NEFF_3_sqrd_phase_cut
        dU_NEFF_TOT_phase_cut=np.sqrt(dU_NEFF_TOT_sqrd_phase_cut)
        
        
        
        
        U_NEFF_TOT_phase_cut=U_NEFF_1_phase_cut+U_NEFF_2_phase_cut+U_NEFF_3_phase_cut


        q_phase_cut=Q_NEFF_TOT_phase_cut/I_NEFF_TOT_phase_cut
        u_phase_cut=U_NEFF_TOT_phase_cut/I_NEFF_TOT_phase_cut

        dq_phase_cut=dQ_NEFF_TOT_phase_cut/I_NEFF_TOT_phase_cut
        du_phase_cut=dU_NEFF_TOT_phase_cut/I_NEFF_TOT_phase_cut

        PD_phase_cut=np.sqrt( (q_phase_cut**2) + (u_phase_cut**2) )
        PD_array.append(PD_phase_cut)
        
        dPD_phase_cut=np.sqrt( (q_phase_cut*dq_phase_cut)**2 + (u_phase_cut*du_phase_cut)**2 )/PD_phase_cut
        dPD_array.append(dPD_phase_cut)


        PA_phase_cut=0.5*math.atan2(u_phase_cut,q_phase_cut)
        PA_array.append(PA_phase_cut)
        
        
        dPA_phase_cut=0.5*(np.sqrt( (u_phase_cut*dq_phase_cut)**2 + (q_phase_cut*du_phase_cut)**2)/ (PD_phase_cut**2))
        dPA_array.append(dPA_phase_cut)
        
        
       
        
        
        
        
        
    PA_array=[np.degrees(i) for i in PA_array]
    dPA_array=[np.degrees(i) for i in dPA_array]
        
    PD_array=[100*i for i in PD_array]
    dPD_array=[100*i for i in dPD_array]

    mean_PD=np.mean(PD_array)
    mean_dPD=np.mean(dPD_array)
    mean_PA=np.mean(PA_array)
    mean_dPA=np.mean(dPA_array)
    print('Mean PD',mean_PD)
    print('Mean PA',mean_PA)
    
    plt.figure()
    plt.title('I vs phase')
    plt.plot(phase_plot,I_array,'.')
    plt.show()
    
    plt.figure()
    plt.title('PA vs phase')
    plt.plot(phase_plot,PA_array,'.')
    plt.xlabel('Phase')
    plt.ylabel('PA (degrees)')
    plt.errorbar(phase_plot,PA_array,yerr=dPA_array,ls='none')
    plt.show()
    
    plt.figure()
    plt.title('PD vs Phase')
    plt.plot(phase_plot,PD_array,'.',label='PD (percent)')
    plt.xlabel('Phase')
    plt.ylabel('PD (percent)')
    plt.errorbar(phase_plot,PD_array,yerr=dPD_array,ls='none')
    plt.legend()
    plt.show()

    results=np.array(tuple(zip(phase_plot, av_phase_err, PD_array, dPD_array, PA_array, dPA_array,I_array)))
    np.savetxt(output_file,results)
    
    return phase_plot, av_phase_err, PD_array, dPD_array, PA_array, dPA_array,I_array
        

