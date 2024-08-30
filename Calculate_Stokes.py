import math
import numpy as np
from astropy.io import fits

def NEFF_STOKES(ixpe_event_file_DU1,
           ixpe_event_file_DU2,
           ixpe_event_file_DU3,ixpe_master_GTI,
           response_dir):
    
    # Assign IXPE min and max energy channels:
    Pmin=51
    Pmax=200
    GTI_master=list(np.loadtxt(str(ixpe_master_GTI)))
            
        # Event file DU1
        
    with fits.open(ixpe_event_file_DU1) as hdu:
        EVENTS_DU1=hdu[1].data
        EVENTS_DU1=EVENTS_DU1[(Pmin<=EVENTS_DU1['PI']) & (EVENTS_DU1['PI']<= Pmax)]
        
        W_MOM_DU1=list(EVENTS_DU1['W_MOM'])
        PI_DU1=list(EVENTS_DU1['PI'])
        q_1=list(EVENTS_DU1['Q'])
        u_1=list(EVENTS_DU1['U'])
        
     

        lc_DU1 =Lightcurve(EVENTS_DU1['TIME'], [1]*len(EVENTS_DU1['TIME']) ,gti=GTI_master)
        lc_DU1.apply_gtis()
        
        
      
        filter_set_DU1=set(lc_DU1.time)
    
    
        
        events_1=[(index,value) for index,value in enumerate(EVENTS_DU1['TIME']) if value in filter_set_DU1]
        indicies_DU1=list([index for index, value in events_1])
        TIMES_DU1=[value for index,value in events_1]
        
        W_MOM_DU1=np.array([W_MOM_DU1[i] for i in indicies_DU1])
        PI_DU1=np.array([PI_DU1[i] for i in indicies_DU1])
        q_1=np.array([q_1[i] for i in indicies_DU1])
        u_1=np.array([u_1[i] for i in indicies_DU1])
        
        
       
        #filtered_another_list = [another_list[index] for index in filtered_indices]

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
        
        W_MOM_DU3=np.array([W_MOM_DU3[i] for i in indicies_DU3])
        PI_DU3=np.array([PI_DU3[i] for i in indicies_DU3])
        q_3=np.array([q_3[i] for i in indicies_DU3])
        u_3=np.array([u_3[i] for i in indicies_DU3])
        
            
  
        
        
        
        #DU1 RESPONSE FILES
    
    
    
    #MRF
    with fits.open(response_dir+'mrf/ixpe_d1_20170101_alpha075_05.mrf') as hdu:
        SPECRESP_MRF_DU1=hdu[1].data
        Aeff_mu_E_DU1=SPECRESP_MRF_DU1['SPECRESP']
        #print(Aeff_mu_E_DU1[0][0])
        
    #ARF
    with fits.open(response_dir+'arf/ixpe_d1_20170101_alpha075_05.arf') as hdu:
        SPECRESP_ARF_DU1=hdu[1].data
        Aeff_E_DU1=SPECRESP_ARF_DU1['SPECRESP']
        
        ENERGY_LO=SPECRESP_MRF_DU1['ENERG_LO']
        ENERGY_HI=SPECRESP_MRF_DU1['ENERG_HI']
        
       
    
        #DU2 RESPONSE FILES
    
    #MRF
    with fits.open(response_dir+'mrf/ixpe_d2_20170101_alpha075_05.mrf') as hdu:
        SPECRESP_MRF_DU2=hdu[1].data
        Aeff_mu_E_DU2=SPECRESP_MRF_DU2['SPECRESP']
        
        
    #ARF
    with fits.open(response_dir+'arf/ixpe_d2_20170101_alpha075_05.arf') as hdu:
        SPECRESP_ARF_DU2=hdu[1].data
        Aeff_E_DU2=SPECRESP_ARF_DU2['SPECRESP']
        
        
        
            #DU3 RESPONSE FILES
    #MRF
    with fits.open(response_dir+'mrf/ixpe_d3_20170101_alpha075_05.mrf') as hdu:
        SPECRESP_MRF_DU3=hdu[1].data
        Aeff_mu_E_DU3=SPECRESP_MRF_DU3['SPECRESP']
    #ARF
    with fits.open(response_dir+'arf/ixpe_d3_20170101_alpha075_05.arf') as hdu:
        SPECRESP_ARF_DU3=hdu[1].data
        Aeff_E_DU3=SPECRESP_ARF_DU3['SPECRESP']
        
        

    #ENERGY_BINS=np.array( list( enumerate( tuple( zip(ENERGY_LO,ENERGY_HI) ),1 ))) 
       
        
    # Assigning response variables to each event   
    Aeff_mu_event_DU1=[Aeff_mu_E_DU1[i-1] for i in PI_DU1]   
    Aeff_mu_event_DU2=[Aeff_mu_E_DU2[i-1] for i in PI_DU2]
    Aeff_mu_event_DU3=[Aeff_mu_E_DU3[i-1] for i in PI_DU3]
    
    Aeff_event_DU1=[Aeff_E_DU1[i-1] for i in PI_DU1]
    Aeff_event_DU2=[Aeff_E_DU2[i-1] for i in PI_DU2]
    Aeff_event_DU3=[Aeff_E_DU3[i-1] for i in PI_DU3]
    
    
    #Calculating weighted Stokes
        #DU1
    I_NEFF_1=np.sum(W_MOM_DU1/Aeff_event_DU1)
    dI_NEFF_1_sqrd=np.sum( (W_MOM_DU1/Aeff_event_DU1)**2 ) 
        #DU2
    I_NEFF_2=np.sum(W_MOM_DU2/Aeff_event_DU2)
    dI_NEFF_2_sqrd=np.sum( (W_MOM_DU2/Aeff_event_DU2)**2 ) 
        #DU3
    I_NEFF_3=np.sum(W_MOM_DU3/Aeff_event_DU3)
    dI_NEFF_3_sqrd=np.sum( (W_MOM_DU3/Aeff_event_DU3)**2 ) 
        #TOT
    dI_NEFF_TOT_sqrd=dI_NEFF_1_sqrd+dI_NEFF_2_sqrd+dI_NEFF_3_sqrd
    dI_NEFF_TOT=np.sqrt(dI_NEFF_TOT_sqrd)
    I_NEFF_TOT=I_NEFF_1+I_NEFF_2+I_NEFF_3
    
    
    
    # Stokes Q
    
        #DU1
    Q_NEFF_1=np.sum( (W_MOM_DU1*q_1) / Aeff_mu_event_DU1)
    dQ_NEFF_1_sqrd=np.sum( ( (W_MOM_DU1*q_1)/Aeff_event_DU1 )**2 )
        #DU2
    Q_NEFF_2=np.sum( (W_MOM_DU2*q_2) / Aeff_mu_event_DU2)
    dQ_NEFF_2_sqrd=np.sum( ( (W_MOM_DU2*q_2)/Aeff_event_DU2 )**2 )
        #DU3
    Q_NEFF_3=np.sum( (W_MOM_DU3*q_3) / Aeff_mu_event_DU3)
    dQ_NEFF_3_sqrd=np.sum( ( (W_MOM_DU3*q_3)/Aeff_event_DU3 )**2 )
        #TOT
    dQ_NEFF_TOT_sqrd=dQ_NEFF_1_sqrd+dQ_NEFF_2_sqrd+dQ_NEFF_3_sqrd
    dQ_NEFF_TOT=np.sqrt(dQ_NEFF_TOT_sqrd)
    Q_NEFF_TOT=Q_NEFF_1+Q_NEFF_2+Q_NEFF_3
    
    
    
   # Stokes U 
        #DU1
    U_NEFF_1=np.sum( (W_MOM_DU1*u_1) / Aeff_mu_event_DU1)
    dU_NEFF_1_sqrd=np.sum( ( (W_MOM_DU1*u_1)/Aeff_event_DU1 )**2 )
        #DU2
    U_NEFF_2=np.sum( (W_MOM_DU2*u_2) / Aeff_mu_event_DU2)
    dU_NEFF_2_sqrd=np.sum( ( (W_MOM_DU2*u_2)/Aeff_event_DU2 )**2 )
        #DU3
    U_NEFF_3=np.sum( (W_MOM_DU3*u_3) / Aeff_mu_event_DU3)
    dU_NEFF_3_sqrd=np.sum( ( (W_MOM_DU3*u_3)/Aeff_event_DU3 )**2 )
        #TOT
    dU_NEFF_TOT_sqrd=dU_NEFF_1_sqrd+dU_NEFF_2_sqrd+dU_NEFF_3_sqrd
    dU_NEFF_TOT=np.sqrt(dU_NEFF_TOT_sqrd)
    U_NEFF_TOT=U_NEFF_1+U_NEFF_2+U_NEFF_3
    
    
    
    
    #print('I_NEFF : ',I_NEFF_TOT)
    #print('dI_NEFF : ',dI_NEFF_TOT)
    #print('Q_NEFF: ',Q_NEFF_TOT)
    #print('dQ_NEFF : ',dQ_NEFF_TOT)
    #print('U_NEFF: ',U_NEFF_TOT)
    #print('dU_NEFF: ',dU_NEFF_TOT)
    
    q=Q_NEFF_TOT/I_NEFF_TOT
    u=U_NEFF_TOT/I_NEFF_TOT
    
    dq=dQ_NEFF_TOT/I_NEFF_TOT
    du=dU_NEFF_TOT/I_NEFF_TOT
    
    PD=100*np.sqrt(q**2+u**2)
    dPD=100*np.sqrt( (q*dq)**2 + (u*du)**2 )/PD
    
    
    PA=np.degrees(0.5*math.atan2(u,q))
    dPA=np.degrees(np.sqrt( (u*dq)**2 + (q*du)**2)/ (2*PD**2))
    
    
    return PD, dPD, PA, dPA
    
    
    
    
    