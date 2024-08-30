from astropy.io import fits
import numpy as np
import math

#Creating a 'master'  GTI file from all DUs of IXPE for a given observation

def master_GTI(data1,data2,data3,output_file):
    
    with fits.open(str(data1)) as hdu1:
        GTI1=list(hdu1[2].data)
        gtistart1=[item[0] for item in GTI1]
        gtiend1=[item[1] for item in GTI1]
        EVENTS_header1=hdu1[1].header
        TSTOP1=EVENTS_header1['TSTOP']
        ngti1=len(GTI1)
    
    
    with fits.open(str(data2)) as hdu2:
        GTI2=list(hdu2[2].data)
        gtistart2=[item[0] for item in GTI2]
        gtiend2=[item[1] for item in GTI2]
        ngti2=len(GTI2)

    with fits.open(str(data3)) as hdu3:
        GTI3=list(hdu3[2].data)
        gtistart3=[item[0] for item in GTI3]
        gtiend3=[item[1] for item in GTI3]
        ngti3=len(GTI3)


    
    
    
    
    beep=0
    counter = 0
    i1 = 0 
    i2 = 0
    i3 = 0
    gtistart=[]
    gtiend=[]
    
    while TSTOP1 - beep > 1e-6:  #while the beep indexer isnt the end of the observation:
        
        gtistart.append(max(gtistart1[i1], gtistart2[i2], gtistart3[i3])) #defining start of current master GTI
        print('gti start {}'.format(gtistart[counter]))
       # print('gti end {}'.format(gtiend))
        gtiend.append(0)
        #print('gti end {}'.format(gtiend))
        gtiend[counter]=TSTOP1 #making the current index of gtiend equal to the end of obs (broken: 
      #  print('gti end counter {}'.format(gtiend[counter]))
        
        
        #first iteration: start is the latest start, stop is the final stop
        
        
        #choosing our end gti
        #if the current du gti ends after current master start (which is first the max start time) then the current master end is the earliest start between that current du end or the current master end
        
        
        if( gtiend1[i1] > gtistart[counter]):
            gtiend[counter]=( min( gtiend1[i1] , gtiend[counter] ))
        else:
            None
        if( gtiend2[i2] > gtistart[counter] ):
            gtiend[counter]=( min( gtiend2[i2] , gtiend[counter] ))
        else:
            None
        if( gtiend3[i3] > gtistart[counter] ):
            gtiend[counter]=( min( gtiend3[i3] , gtiend[counter] ))
        else:
            None
    
        
       
        #moving to the next current du gti if it has an end before the master end
        
        
        if gtiend1[i1] - gtiend[counter] < 1e-6 and i1<ngti1:
            i1 = i1 + 1
        else:
            None
        if gtiend2[i2] - gtiend[counter] < 1e-6 and i2<ngti2:
            i2 = i2 + 1
        else:
            None
        if gtiend3[i3] - gtiend[counter] < 1e-6 and i3<ngti3:
            i3 = i3 + 1
        else:
            None
            
            
        print(i1)
        print(i2)
        print(i3)
            
            
    #moving to the next current du gti if it has an end time before the latest current start
        if i1<ngti1 and i2<ngti2 and i3<ngti3:
            tlatest=max(gtistart1[i1],gtistart2[i2],gtistart3[i3])
        
      #  else:
      #      None 
            if gtiend1[i1] - tlatest < 1e-6 and i1<ngti1:
                i1 = i1 + 1
            else:
                None
            if gtiend2[i2] - tlatest < 1e-6 and i2<ngti2:
                i2 = i2 + 1
            else:
                None
            if gtiend3[i3] - tlatest < 1e-6 and i3<ngti3:
                i3 = i3 + 1
            else:
                None
                
        else:
            None

            
        #print(i1)
       # print(gtistart1[i1])
        beep=gtiend[counter]   
        counter=counter+1
        
        
    #return [gtistart,gtiend]

    master_gti=[list(x) for x in zip(gtistart,gtiend)]
    np.savetxt(str(output_file),master_gti)
    #print((master_gti))
    #print(len(master_gti))
    
