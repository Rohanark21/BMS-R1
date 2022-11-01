import numpy as np
import pandas as pd
from panda_init import initSPKF
from panda_iter import iterSPKF


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
spkfData=initSPKF(voltage[0],T,SigmaX0,SigmaV,SigmaW,model)

print(spkfData)