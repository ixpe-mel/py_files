def counts_phase(DU123,GTI,Pmin,Pmax,bin_num,v,v_dot,v_double_dot,t_0):
    
    import sys
    import numpy as np
    from stingray import Lightcurve
    sys.path.append('/home/c2032014/PhD/py_files')
    from astropy.io import fits
    import spin_phase_calculator as phase
    
    GTI=list(np.loadtxt(GTI))
    
    with fits.open(DU123) as hdu:
        data=hdu[1].data
        header=hdu[1].header
        data=data[(Pmin<=data['PI']) & (data['PI']<= Pmax)]
        TSTART=header['TSTART']
        TIMES123=data['TIME']
    
    lc_DU123 = Lightcurve(TIMES123, [1]*len(TIMES123) ,gti=GTI)
    lc_DU123.apply_gtis()
    filter_set_DU123=set(lc_DU123.time)
    
    events_1=[(index,value) for index,value in enumerate(TIMES123) if value in filter_set_DU123]
    indicies_DU1=list([index for index, value in events_1])
    TIMES123=[value for index,value in events_1]
    
    phase_per_photon=phase.phase_v_double_dot(v,v_dot,v_double_dot,np.array(TIMES123),t_0)
    
    b=np.histogram(phase_per_photon,bin_num)
    hist_counts=np.array(b[0])
    bin_centers = 0.5 * (b[1][:-1] + b[1][1:])
    bin_length=b[1][1:]-b[1][:-1]
    #plt.figure()
    #lt.hist(phase_per_photon,80)
    #plt.title('Histogram of phase (photons)')
    #plt.xlabel('Phase')
    #plt.ylabel('Counts')
    #plt.plot(bin_centers,b[0],'.')
    return hist_counts , bin_centers , bin_length

    
