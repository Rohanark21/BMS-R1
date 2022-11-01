import numpy as np
import pandas as pd

def SOCfromOCVtemp(ocv,temp,model):  
# This function returns an estimate of soc from a fully rested open-circuit-voltage of an LiPB cell
    
    ocvcol=ocv
    OCV=model['OCV']                     # force ocv to be col-vector
    SOC0=model['SOC0']
    SOCrel=model['SOCrel']
    if np.isscalar(temp): 
        tempcol = temp*np.ones(len(ocvcol))          # replicate temperature for all ocvs
    else:
        tempcol = temp()                             # force to be col vector

    diffOCV = OCV[1] - OCV[0]
    soc = np.zeros((len(ocvcol)))
    I1=ocvcol[ocvcol <= OCV[0]].index
    I2=ocvcol[ocvcol >= OCV[-1]].index
    I3=ocvcol[(ocvcol > OCV[0]) & (ocvcol < OCV[-1])].index
    I6 = np.isnan(ocvcol)
    
    #  for socs lower than lowest voltage, extrapolate off low end of table
    dz = (SOC0[1]+tempcol*SOCrel[1]) - (SOC0[0]+tempcol*SOCrel[0])
    soc[I1] = (ocvcol[I1] - OCV[0])*dz[I1] / diffOCV + SOC0[0] + tempcol[I1]*SOCrel[0]
    
    #  for socs higher than highest voltage, extrapolate off high end of table
    dz = (SOC0[-1]+tempcol*SOCrel[-1]) - (SOC0[-2]+tempcol*SOCrel[-2])
    soc[I2] = (ocvcol[I2]-OCV[-1])*dz[I2]/diffOCV + SOC0[-1]+tempcol[I2]*SOCrel[-1]
    
    #  for normal soc range...manually interpolate (10x faster than "interp1")
    I4=(ocvcol[I3]-OCV[0])/diffOCV
    I5=(np.floor(I4))
    I5 = np.array(list(map(int, I5)))
    soc[I3] = SOC0[I5]*(1 - (I4 - I5)) + SOC0[I5 + 1]*(I4 - I5)
    soc[I3] = soc[I3] + tempcol[I3]*(SOCrel[I5]*(1 - (I4 - I5)) + SOCrel[I5 + 1]*(I4 - I5))
    soc[I6] = 0                                                  # replace NaN OCVs with zero SOC... 03/23/10
    soc = np.reshape(soc,np.size(ocvcol))
    return soc




    