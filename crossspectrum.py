# Make a cross spectrum between two fits files

# Result saves out as [freq,real,im,real*freq,im*freq,cs_err]

def crossspectrum(file1,file2,name,gti,bin_length,seg_length,output_file):
    
    Pmin=51
    Pmax=200
    
    GTI=list(np.loadtxt(str(gti)))
    
    with fits.open(str(file1)) as hdu:
            data_1=hdu[1].data  #reading in DU1+DU2
            TIME1=data_1.field('TIME')
         
    with fits.open(str(file2)) as hdu2:
            data_2=hdu2[1].data #reading in DU3
            TIME3=data_3.field('TIME')

           
    #PI channel/energy index
    index_energy_1=list(locate(data_1.field('PI'), lambda x: Pmin < x < Pmax))  
    data_1=data_1[index_energy_1]
    
    index_energy_2=list(locate(data_2.field('PI'), lambda x: Pmin < x < Pmax))
    data_2=data_2[index_energy_2]

    #Lightcurve

    lightcurve_1=Lightcurve.make_lightcurve(TIME1,dt=bin_length,gti=GTI)
    lightcurve_1.apply_gtis()

    lightcurve_2=Lightcurve.make_lightcurve(TIME2,dt=bin_length,gti=GTI)
    lightcurve_2.apply_gtis()
 

    plt.figure()
    plt.plot(lightcurve_1.time,lightcurve_1.counts)
    plt.title('File 1: Lightcurve')
    plt.show()
    plt.figure()
    plt.plot(lightcurve_3.time,lightcurve_3.counts)
    plt.title('File 2: Lightcurve')
    plt.show()
    
    #Cross spec 
    
    cs = AveragedCrossspectrum.from_lightcurve(lightcurve_12,lightcurve_3,seg_length,norm='frac')
    
    real=cs.power.real
    im=cs.power.imag
    cs_err=cs.power_err

   #Plotting the real part against fourier frequency

    fig, ax1 = plt.subplots(1,1,figsize=(9,6))
        
    ax1.plot(cs.freq, cs.power.real*avg_cs.freq,'.', color='green')
    ax1.errorbar(avg_cs.freq,avg_cs.power.real*avg_cs.freq,yerr=cs_err*avg_cs.freq)
    ax1.set_xlabel("Frequency (Hz)")
    ax1.set_title('Cross spec')
    ax1.set_ylabel("Power x Fourier frequency")
    ax1.set_yscale('log')
    ax1.set_xscale('log')
    ax1.tick_params(axis='x', labelsize=16)
    ax1.tick_params(axis='y', labelsize=16)
    ax1.tick_params(which='major', width=1.5, length=7)
    ax1.tick_params(which='minor', width=1.5, length=4)
    for axis in ['top', 'bottom', 'left', 'right']:
        ax1.spines[axis].set_linewidth(1.5)
        plt.show()

        
    result=np.array(tuple(zip(cs.freq,cs.real,cs.im,cs.real*freq,cs.imag*freq,cs_mean)))
    np.savetxt(output_file,result)