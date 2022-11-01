import numpy as np
import scipy.io
import pandas as pd
from initSPKF_test import initSPKF
from SOCfromOCVtemp_test import SOCfromOCVtemp
from OCVfromSOCtemp_test import OCVfromSOCtemp
import matplotlib.pyplot as plt
np.set_printoptions(formatter={'float_kind':'{:.4f}'.format})

matdata=scipy.io.loadmat('PANdata_P25.mat', simplify_cells=True)
rawdata = pd.DataFrame(matdata['DYNData'])
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
    modeldata = spkfData['model']
  
    #  Load the cell model parameters
    Q  = modeldata['model']['QParam'][5]
    G  = modeldata['model']['GParam'][5]
    M  = modeldata['model']['MParam'][5]
    M0 = modeldata['model']['M0Param'][5]
    RC = np.transpose(np.exp((-deltat)/abs(modeldata['model']['RCParam'][5])))
    R  = np.transpose(modeldata['model']['RParam'][5])
    R0 = modeldata['model']['R0Param'][5]
    eta = modeldata['model']['etaParam']
    if ik<0:
        ik=ik*eta
        ik=ik[5]

    
    #  Get data stored in spkfData structure
    I = spkfData['priorI']
    SigmaX = spkfData['SigmaX']
    xhat = spkfData['xhat']
    Nx = spkfData['Nx']
    Nw = spkfData['Nw']
    Nv = spkfData['Nv']
    Na = spkfData['Na']
    Snoise = spkfData['Snoise']
    Wc = spkfData['Wc']
    irInd = spkfData['irInd']
    hkInd = spkfData['hkInd']
    zkInd = spkfData['zkInd']
    if (abs(ik) > Q/100):
        spkfData['signIk'] = np.sign(ik)
    signIk = spkfData['signIk']
    
    #  Step 1a: State estimate time update
                # - Create xhatminus augmented SigmaX points
                # - Extract xhatminus state SigmaX points
                # - Compute weighted average xhatminus(k)

    #  Step 1a-1: Create augmented SigmaX and xhat
    sigmaXa = np.linalg.cholesky(SigmaX)
    # if p>0:
    #   print('Cholesky error.  Recovering...\n')
    #   theAbsDiag = abs(np.diag(SigmaX))
    #   sigmaXa = np.diag(max(np.sqrt(theAbsDiag),np.sqrt(spkfData['SigmaW']))))
    set1 = np.zeros((3,2))
    set2 = np.zeros((2,3))
    set3 = np.hstack((set2,Snoise))
    sigmaXa = np.hstack((sigmaXa,set1))
    sigmaXa = np.vstack((sigmaXa,set3))
    xhata =np.append(xhat,np.zeros([Nw+Nv,1]))
    #  NOTE: sigmaXa is lower-triangular
    
    #  Step 1a-2: Calculate SigmaX points
    Xa= np.transpose(np.tile(xhata,(11,1)))
    sigmasum=np.zeros((5,1))
    sigmasum=np.hstack((sigmasum,sigmaXa))
    sigmasum=np.hstack((sigmasum,-sigmaXa))
    Xb=sigmasum*spkfData['h']
    Xa=np.add(Xa,Xb)
    
     #  Calculate new states for all of the old state vectors in xold.  
    def stateEqn(xold,current,xnoise,RC,irInd,hkInd,zkInd,G,deltat,Q):
        current = current + xnoise # noise adds to current
        xnew = 0*xold
        xnew[irInd,] = RC*xold[irInd,] + [1-RC]*current
        Ah = np.exp(-abs(current*G*deltat/(3600*Q)))  # hysteresis factor
        xnew[hkInd,] = Ah*xold[hkInd,] + (Ah-1)*np.sign(current)
        xnew[zkInd,] = xold[zkInd,] - current/3600/Q
        xnew[hkInd,] = np.minimum(1,np.maximum(-1,xnew[hkInd,]))
        xnew[zkInd,] = np.minimum(1.05,np.maximum(-0.05,xnew[zkInd,]))
        return xnew

    #  Step 1a-3: Time update from last iteration until now
    #  stateEqn(xold,current,xnoise)
    Xx = stateEqn(Xa[0:3,],I,Xa[3,],RC,irInd,hkInd,zkInd,G,deltat,Q)
    xhat = np.matmul(Xx,spkfData['Wm'])
    
    #  Step 1b: Error covariance time update
    # - Compute weighted covariance sigmaminus(k)
    Xs = Xx - np.transpose(np.tile(xhat,(2*Na+1,1)))
    SigmaX = np.matmul(Xs,np.diag(Wc))
    SigmaX = np.matmul(SigmaX,np.transpose(Xs))
    
     # Calculate cell output voltage for all of state vectors in xhat
    def outputEqn(xhat,current,ynoise,T,modeldata,zkInd,M,hkInd,M0,signIk,R,irInd,R0):
        yhat = OCVfromSOCtemp(xhat[zkInd,],T,modeldata)
        yhat = yhat + M*xhat[hkInd,] + M0*signIk
        yhat = yhat - R*xhat[irInd,] - R0*current + ynoise
        return yhat
    
    #  Step 1c: Output estimate
    # - Compute weighted output estimate yhat(k)
    I = ik 
    yk = vk
    Y = outputEqn(Xx,I,Xa[4,],Tk,modeldata,zkInd,M,hkInd,M0,signIk,R,irInd,R0)
    yhat = np.matmul(Y,spkfData['Wm'])
    
    #  Step 2a: Estimator gain matrix
    Ys = Y - np.transpose(np.tile(yhat,(2*Na+1,1)))
    SigmaXY = np.matmul(Xs,np.diag(Wc))
    SigmaXY = np.matmul(SigmaXY,np.transpose(Ys))
    SigmaY = np.matmul(Ys,np.diag(Wc))
    SigmaY = np.matmul(SigmaY,np.transpose(Ys))
    L = SigmaXY/SigmaY

    #  Step 2b: State estimate measurement update
    r = yk - yhat                                                               # residual.  Use to check for sensor errors...
    if r**2 > 100*SigmaY:
        L[:,0]=0.0
    xhat = xhat + np.squeeze(L)*r
    xhat[zkInd] = np.minimum(1.05,np.maximum(-0.05,xhat[zkInd]))
    xhat[hkInd] = np.minimum(1,np.maximum(-1,xhat[hkInd]))

    #  Step 2c: Error covariance measurement update
    SigmaX = SigmaX - np.matmul(np.matmul(L,SigmaY),np.transpose(L))
    [V,S,U] = np.linalg.svd(SigmaX)
    HH = V*S@np.transpose(V)
    SigmaX = (SigmaX + np.transpose(SigmaX) + HH + np.transpose(HH)) / 4      # Help maintain robustness
    
    #  Q-bump code
    if r**2 > 4*SigmaY:                                                       # bad voltage estimate by 2-SigmaX, bump Q
        print('Bumping sigmax\n')
        SigmaX[zkInd,zkInd] = SigmaX[zkInd,zkInd]*spkfData['Qbump']
        
    #  Save data in spkfData structure for next time...
    spkfData['priorI'] = ik
    spkfData['SigmaX'] = SigmaX
    spkfData['xhat'] = xhat
    zk = xhat[zkInd]
    zkbnd = 3*np.sqrt(SigmaX[zkInd,zkInd])


    # update waitbar periodically, but not too often (slow procedure)
    if k % 5000 == 0:
        print('  Completed {0} out of {1} iterations...\n'.format(k,len(voltage)))
    ind = np.where(abs(soc-sochat)>socbound)
    ind=np.transpose(ind)
    sochat[k]=zk
    socbound[k]=zkbnd
    
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