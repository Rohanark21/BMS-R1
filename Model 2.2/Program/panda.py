import numpy as np
import pandas as pd
from panda_init import initSPKF
from iterSPKF import iterSPKF
import matplotlib.pyplot as plt
np.set_printoptions(formatter={'float_kind':'{:.4f}'.format})

script1=pd.read_excel("PANdata_P25.xlsx")
model=pd.read_excel("PANmodel.xlsx")

T=25
eta=script1['eta']
voltage=script1['voltage']
current=script1['current']
time=script1['time']
soc=script1['soc']
deltat = time[1]-time[0]
time    = time-time[0]  

# Reserve storage for computed results, for plotting
sochat = np.zeros(np.size(soc))
socbound = np.zeros(np.size(soc))

# Covariance values
SigmaX0 = np.diag(np.array([100.0,0.01,0.001]))              # uncertainty of initial state
SigmaV = 0.3                                                 # Uncertainty of voltage sensor, output equation
SigmaW = 4.0                                                 # Uncertainty of current sensor, state equation

# Create spkfData structure and initialize variables using first voltage measurement and first temperature measurement
spkfData=initSPKF(voltage[0],T,SigmaX0,SigmaV,SigmaW,script1)

print('Please be patient. This code will take several hours to execute.\n' % ())
for k in range(0,len(voltage)):                              #len(voltage)
    vk = voltage[k]
    ik = current[k]
    Tk = T
    
    # Update SOC (and other model states)
    sochat[k],socbound[k],spkfData = iterSPKF(vk,ik,Tk,deltat,spkfData)
    # update waitbar periodically, but not too often (slow procedure)
    if k % 5000 == 0:
        print('  Completed {0} out of {1} iterations...\n'.format(k,len(voltage)))
    ind = np.where(abs(soc-sochat)>socbound)
    ind=np.transpose(ind)
    
print('RMS SOC estimation error = %g%%' % (np.sqrt(np.mean((100 * (soc - sochat)) ** 2))))
print('And Percent of time error outside bounds = %g%%\n' % (len(ind)/len(soc)*100))