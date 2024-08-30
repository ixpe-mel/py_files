from astropy.table import Table
from astropy.io import fits
import numpy as np

# <first ixpe fits file> <second ixpe fits file> <output file name>

# Stiches together IXPE files: Only use GTI concatenation if GTIs come from same DU

def concatenate_ixpe_fits(file1, file2, output_file):
    with fits.open(file1) as hdu:
        data1=Table(hdu[1].data)
        GTI1=Table(hdu[2].data) 
        prihdr = hdu[0].header # Define primary hdu header
        events_header=hdu[1].header
        GTI_header=hdu[2].header
        prihdu = fits.PrimaryHDU(header=prihdr) #Define primary hdu
        #print(type(GTI1))
        #print(GTI1)
        
    with fits.open(file2) as hdu2:
        data2=Table(hdu2[1].data)
        GTI2=Table(hdu2[2].data)
    
    # Concatenate the tables containing the hdu data

    concatenated_events_table=np.hstack([data1,data2])  
    concatenated_GTI_table=np.hstack([GTI1,GTI2])
    
    events_hdu=fits.BinTableHDU.from_columns(concatenated_events_table,name='EVENTS',header=events_header)
    GTI_hdu=fits.BinTableHDU.from_columns(concatenated_GTI_table,name='GTI',header=GTI_header)
    
    thdulist = fits.HDUList([prihdu, events_hdu,GTI_hdu]) #Defining the list of hdus
    thdulist.writeto(output_file, overwrite=True)
