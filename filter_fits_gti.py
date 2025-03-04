#We usually need to make a fourier product or lightcurve in order to 
#filter of GTIs, so here we can filter a dataset over GTIs using stingray

import numpy as np
from astropy.io import fits
from stingray import Lightcurve

#We want the events 'data' to be already loaded in and cleaned

def filter_fits_gti(data,ixpe_master_GTI):

    GTI_master=list(np.loadtxt(str(ixpe_master_GTI)))
    W_MOM=data['W_MOM']
    PI=data['PI']
    q=data['Q']
    u=data['U']
    PHI=data['PHI']


    #Using stingray to filter events over the GTI
    lc = Lightcurve(data['TIME'], [1]*len(data['TIME']) ,gti=GTI_master)
    lc.apply_gtis()
    filter_set=set(lc.time)

    events=[(index,value) for index,value in enumerate(data['TIME']) if value in filter_set]
    #print(events_1[:100])
    indicies=list([index for index, value in events])
    TIMES=[value for index,value in events]
    
    W_MOM=np.array([W_MOM[i] for i in indicies])
    PI=np.array([PI[i] for i in indicies])
    q=np.array([q[i] for i in indicies])
    u=np.array([u[i] for i in indicies])
    PHI=np.array([PHI[i] for i in indicies])

    
    return np.array(TIMES),np.array(W_MOM),np.array(PI),np.array(q),np.array(u),np.array(PHI)


    #W_MOM_DU1=np.array([W_MOM_DU1[i] for i in indicies_DU1])
    #PI_DU1=np.array([PI_DU1[i] for i in indicies_DU1])
    #q_1=np.array([q_1[i] for i in indicies_DU1])
    #u_1=np.array([u_1[i] for i in indicies_DU1])