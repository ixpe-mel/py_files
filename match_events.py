#takes in an ixpe file and a corrected binarty orbit file and adds a column of the uncorrected times to the corrected file, and saves it to an output file

def match_events(input_uncorrected_file,input_corrected_file,output_correction_added_file):
    
    with fits.open(input_corrected_file) as hdu:
        data_corr=hdu[1].data
        TIMES_corr=data_corr['TIME']
        PI_corr=data_corr['PI']

    
    
    with fits.open(input_uncorrected_file) as hdu:
        data=hdu[1].data
        TIMES=data['TIME']
        PI=data['PI']
    
    
    
    uncor_matched_times=[]   
    index=0
    for i in range(len(PI_corr)):
 
        condition_met=False
        while not condition_met:
            current_cor=cor_PI[i]
            #print('corrected PI',current_cor)
            current_uncor=uncor_PI[index]
            #print('current_uncor',current_uncor)
            
            if current_cor==current_uncor:
                condition_met=True
                #print('match')
                uncor_matched_times.append(uncor_times[index])
            index+=1
    #return uncor_matched_times
    
    # Step 1: Load the FITS file
    file_path = input_corrected_file
    hdul = fits.open(file_path)
    data = Table(hdul[1].data)
    
    prihdr = hdu[0].header
    prihdu = fits.PrimaryHDU(header=prihdr) #Define primary hdu
    
    data.rename_column('TIME','TIME_corr')
    events_header=hdul[1].header



    # Step 3: Create a new column with the NumPy array data
    new_column_data = np.array(uncor_matched_times) # Replace with your actual data
    new_column = Column(new_column_data, name='TIME')

# Step 4: Add the new column to the table
    data.add_column(new_column)

    cols=data.colnames

    col1_index=cols.index('TIME_corr')
    col2_index=cols.index('TIME')


    cols[col1_index], cols[col2_index] = cols[col2_index], cols[col1_index]
    data = data[cols]

    new_hdu = fits.BinTableHDU(data,header=events_header,name='EVENTS')

# Step 6: Write the updated HDU list to a new FITS file
    new_hdul = fits.HDUList([prihdu, new_hdu])  # Include the primary HDU and the new table HDU
    new_hdul.writeto('/home/c2032014/PhD/binarycorr/adam/02004001/event_l2/ixpe02004001_det1_evt2_timecorrfiltered_plus_uncorr.fits', overwrite=True)
    
    return uncor_matched_times


            