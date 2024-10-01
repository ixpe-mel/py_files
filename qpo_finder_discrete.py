#Discretised version of QPO finder
import sys
sys.path.append('/home/c2032014/PhD/py_files/')
import load_and_clean as lac
import py_files.rms_phase_calculator_NEW as rpc
import rms_normalisation as rn
import numpy as np

def qpo_finder(file1,file2,output_file,gti,bin_length,seg_length,fmin,fmax,Pmin,Pmax,mod_bin_number):

    #load data and clean
    data_1,header_1=lac.load_and_clean(file1,Pmin,Pmax)
    data_2,header_2=lac.load_and_clean(file2,Pmin,Pmax)

    #calculate rms normalisation
    norm_factor,lc_ref=rn.rms_normalisation(data_1,data_2,bin_length,seg_length,fmin,fmax,gti)

    #Calculate G, subtracting spurious polarisation

    frac_rms,phase=rpc.rms_phase_calc(data_1,lc_ref,bin_length,seg_length,fmin,fmax,norm_factor,mod_bin_number)
    return frac_rms,phase