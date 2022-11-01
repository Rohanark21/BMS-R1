import numpy as np
import scipy.io
from initSPKF import initSPKF
from iterSPKF import iterSPKF
import matplotlib.pyplot as plt
np.set_printoptions(formatter={'float_kind':'{:.4f}'.format})

matdata=scipy.io.loadmat('PANdata_P25.mat', simplify_cells=True)
rawdata = matdata['DYNData']
modeldata=scipy.io.loadmat('PANmodel.mat', simplify_cells=True)

T=25
eta = rawdata['eta']
time = rawdata['script1']['time']
deltat = time[1]-time[0]
time    = time-time[0]                                      # start time at 0
current = rawdata['script1']['current']                     # discharge > 0; charge < 0.
voltage = rawdata['script1']['voltage']
soc     = rawdata['script1']['soc']

# Reserve storage for computed results, for plotting
sochat = np.zeros(np.size(soc))
socbound = np.zeros(np.size(soc))

# Covariance values
SigmaX0 = np.diag(np.array([100.0,0.01,0.001]))              # uncertainty of initial state
SigmaV = 0.3                                                 # Uncertainty of voltage sensor, output equation
SigmaW = 4.0                                                 # Uncertainty of current sensor, state equation

# Create spkfData structure and initialize variables using first voltage measurement and first temperature measurement
spkfData=initSPKF(voltage[0],T,SigmaX0,SigmaV,SigmaW,modeldata)

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

### plotting program ###
plt.subplot(1,2,1)
plt.plot(time / 60,100 * sochat,time / 60,100 * soc,linewidth=0.5)
plt.plot(time/60,100*(sochat+socbound),time/60,100*(sochat-socbound),color='green',linewidth=0.5)
plt.title('SOC estimation using SPKF')
plt.grid('on')
plt.xlabel('Time (min)')
plt.ylabel('SOC (%)')
plt.legend(('Estimate','Truth','Bounds'))

plt.subplot(1,2,2)
plt.plot(time / 60,100 * (soc - sochat),linewidth=0.5)
plt.plot(time/60,100*socbound,time/60,-100*socbound,color='m',linewidth=0.5)
plt.title('SOC estimation errors using SPKF')
plt.xlabel('Time (min)')
plt.ylabel('SOC error (%)')
plt.ylim(np.array([- 4,4]))
plt.legend(('Estimation error','Bounds'))
plt.grid('on')
plt.show()
print("zhala")