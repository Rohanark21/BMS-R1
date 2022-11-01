import numpy as np
from SOCfromOCVtemp import SOCfromOCVtemp
from OCVfromSOCtemp import OCVfromSOCtemp


def iterSPKF(vk,ik,Tk,deltat,spkfData):
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
  Xa= np.tile(xhata,(11,1))
  Xa = Xa.transpose()
  sigmasum=np.zeros((5,1))
  sigmasum=np.hstack((sigmasum,sigmaXa))
  sigmasum=np.hstack((sigmasum,-sigmaXa))
  Xb=sigmasum*spkfData['h']
  Xa=np.add(Xa,Xb)
 
  
  #  Step 1a-3: Time update from last iteration until now
  #  stateEqn(xold,current,xnoise)
  Xx = stateEqn(Xa[0:3,],I,Xa[3,],RC,irInd,hkInd,zkInd,G,deltat,Q)
  xhat = np.matmul(Xx,spkfData['Wm'])
  
  #  Step 1b: Error covariance time update
  # - Compute weighted covariance sigmaminus(k)
  Xs = Xx - np.transpose(np.tile(xhat,(2*Na+1,1)))
  SigmaX = np.matmul(Xs,np.diag(Wc))
  SigmaX = np.matmul(SigmaX,np.transpose(Xs))
  
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
  return zk,zkbnd,spkfData

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

  # Calculate cell output voltage for all of state vectors in xhat
def outputEqn(xhat,current,ynoise,T,modeldata,zkInd,M,hkInd,M0,signIk,R,irInd,R0):
  yhat = OCVfromSOCtemp(xhat[zkInd,],T,modeldata)
  yhat = yhat + M*xhat[hkInd,] + M0*signIk
  yhat = yhat - R*xhat[irInd,] - R0*current + ynoise
  return yhat


