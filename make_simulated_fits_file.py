# Writing a code that adds the simulated modulation angle column to a
#  fits file to replace the real data

import numpy as np
import sys
from astropy.table import Table, Column
sys.path.append('/home/c2032014/py_files')
import load_and_clean as lac


def create_sim_fits_file(file,Pmin,Pmax,simulated_mod_angles,output_file):
    
    #load true fits file
    data,header=lac.load_and_clean(file,Pmin,Pmax)
    hdutable=Table(data)
    print(hdutable['PHI'])
    #load simulated modulation angles

    sim_mod_angles=np.loadtxt(simulated_mod_angles)
    # Add the simulated modulation angles to the data
    
    # Change the name of the column in the FITS file
    #data.rename_column('phi', 'phi_true')
    #col_name = 'phi'  
    hdutable['PHI'].name = 'PHI_TRUE'

    # Add the simulated modulation angles to the data
    new_column=Column(sim_mod_angles,name='PHI')
    hdutable['PHI']=new_column

    # Create a new FITS file with the updated data
    hdu = fits.BinTableHDU(hdutable, header=header)
    hdu.writeto(output_file, overwrite=True)