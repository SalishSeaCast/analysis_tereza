import numpy as np

def atm_co2_gauss_trunc(year, yearday):
    
    y2 = 0
    x = yearday 

    LR_slope   =  2.149     
    LR_int     =  -3929.359 
    ctr        =  161.898  
    amp        =  7.083     
    wid        =  44.703    
    ctr2       =  218.832    
    amp2       =  -19.004   
    wid2       =  87.8836   
    ctr3       =  199.430   
    amp3       =  8.026     
    wid3       =  -185.920 
    
    y2 = y2 + amp * np.exp( -((x - ctr)/wid)**2)\
    + amp2 * np.exp( -((x - ctr2)/wid2)**2)\
    + amp3 * np.exp( -((x - ctr3)/wid3)**2)
    
    finval = (year+yearday/365)*LR_slope+LR_int + y2
    
    return finval