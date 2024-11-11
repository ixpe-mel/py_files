#Stores all the potential models that can be used to fit the rms/phase
import numpy as np

def line(x,c):
        y=(x*0)+c  
        return y
        
def sin_sum(x,A_01,A1,delta1,A2,delta2):
        y=A_01 + A1*(np.sin(4*np.array(x)+delta1)) +A2*(np.sin(2*np.array(x)+delta2))
        return y
    
    
def sin_90(x,A_0,A,delta):
        y=A_0 + A*(np.sin((4*(np.array(x)+delta))))
        return y
    
def sin_180(x,A_0,A,delta):
        y=A_0 + A*(np.sin((2*(np.array(x)+delta))))
        return y