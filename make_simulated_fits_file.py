# Writing a code that adds the simulated modulation angle column to a
#  fits file to replace the real data

import numpy as np
import sys
from astropy.io import fits
from astropy.table import Table, Column
sys.path.append('/home/c2032014/py_files')
import load_and_clean as lac


def create_sim_fits_file(file,Pmin,Pmax,simulated_mod_angles,output_file):
    
    #load true fits file
    with fits.open(file) as hdu:
        
        
        data,events_header,GTI,GTI_header,prihdu=lac.load_and_clean(file,51,200)
        data=Table(data)
        GTI=Table(GTI) 
        #prihdr = hdu[0].header # Define primary hdu header
        #events_header=hdu[1].header
        #GTI_header=hdu[2].header
        #prihdu = fits.PrimaryHDU(header=prihdr) #Define primary hdu
    
    
    #print(hdutable['PHI'])
    #load simulated modulation angles

    sim_mod_angles=np.loadtxt(simulated_mod_angles)
    # Add the simulated modulation angles to the data
    
    # Change the name of the column in the FITS file
    #data.rename_column('phi', 'phi_true')
    #col_name = 'phi'  
    data['PHI'].name = 'PHI_TRUE'

    # Add the simulated modulation angles to the data
    new_column=Column(sim_mod_angles,name='PHI')
    data['PHI']=new_column

    # Create a new FITS file with the updated data
    events_hdu = fits.BinTableHDU(data,name='EVENTS', header=events_header)
    GTI_hdu=fits.BinTableHDU(GTI,name='GTI',header=GTI_header)
    
    thdulist = fits.HDUList([prihdu, events_hdu,GTI_hdu])

    thdulist.writeto(output_file, overwrite=True)


def add_fits_column(fits_file,text_file,column_name,output_file):
    data,events_header,GTI,GTI_header,prihdu=lac.load_and_clean(fits_file,51,200)
    data=Table(data)
    new_column=Column(np.loadtxt(text_file),name=column_name)
    data[column_name]=new_column
        #print(data)
        #print(data[column_name])
        #hdu[1].data=data
    events_hdu = fits.BinTableHDU(data,name='EVENTS', header=events_header)
    GTI_hdu=fits.BinTableHDU(GTI,name='GTI',header=GTI_header)
    thdulist = fits.HDUList([prihdu, events_hdu,GTI_hdu])
    thdulist.writeto(output_file,overwrite=True)
    