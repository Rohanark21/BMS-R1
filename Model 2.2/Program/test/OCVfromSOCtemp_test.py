import numpy as np

def OCVfromSOCtemp(soc,temp,modeldata):
#  This function returns the fully rested open-circuit-voltage of an LiPB cell given its soc.
    
    soccol=np.squeeze(np.array([soc])) 
    SOC = modeldata['model']['SOC']                                     # force to be col vector... 03/24/10
    OCV0 = modeldata['model']['OCV0']                                   # force to be col vector... 03/24/10
    OCVrel = modeldata['model']['OCVrel']                               # force to be col vector... 03/24/10
    if np.isscalar(temp): 
        tempcol = temp*np.ones(len(soccol))                             # replicate for all socs
    else:
        tempcol = temp()                                                # force to be col vector
    
    diffSOC=SOC[1]-SOC[0]
    ocv=np.zeros(len(soccol))
    I1=np.where(soccol <= SOC[1])                                        # if ~isempty(I1), disp('low soc'); end
    I2=np.where(soccol >= SOC[-1])                                       # if ~isempty(I2), disp('high soc'); end
    I3=np.where((soccol > SOC[0]) & (soccol < SOC[-1]))
    I6=np.isnan(soccol)
    # print("I3")
    # print(I3)
    # print(np.shape(I3))
    #  for voltages less than 0% soc... 07/26/06
    #  extrapolate off low end of table (for SOC(1) < 0... 03/23/10)
    if len(I1) !=0:
        dv = (OCV0[1]+tempcol*OCVrel[1]) - (OCV0[0]+tempcol*OCVrel[0])
        ocv[I1]= (soccol[I1]-SOC[0])*dv[I1]/diffSOC + OCV0[0]+tempcol[I1]*OCVrel[0]   
    
    #  for voltages greater than 100% soc... 07/26/06
    #  extrapolate off high end of table (for SOC(end) > 1... 03/23/10)
    if len(I2) !=0:
        dv = (OCV0[-1]+tempcol*OCVrel[-1]) - (OCV0[-2]+tempcol*OCVrel[-2])
        ocv[I2] = (soccol[I2]-SOC[-1])*dv[I2]/diffSOC + OCV0[-1]+tempcol[I2]*OCVrel[-1]
    
    #  for normal soc range...
    #  manually interpolate (10x faster than "interp1")
    I4=(soccol[I3]-SOC[0])/diffSOC                                            # for SOC(1) < 0... 03/23/10
    I5=(np.floor(I4))
    I5 = np.array(list(map(int,I5)))
    I45=I4-I5
    omI45=1-I45
    ocv[I3] = OCV0[I5]*omI45 + OCV0[I5+1]*I45
    ocv[I3] = ocv[I3] + tempcol[I3]*(OCVrel[I5]*omI45 + OCVrel[I5+1]*I45)
    ocv[I6] = 0                                                            # replace NaN SOCs with zero voltage... 03/23/10
    ocv = np.reshape(ocv,np.size(soc))
    return ocv

