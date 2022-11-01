import numpy as np
from SOCfromOCVtemp_test import SOCfromOCVtemp

def initSPKF(voltage,temperature,SigmaX0,SigmaV,SigmaW,model):

#   Initial state description
  spkfData = dict()     
  ir0 = 0                             
  hk0 = 0                           
  SOC0 = SOCfromOCVtemp(voltage,temperature,model)
  spkfData['irInd'] = 0
  spkfData['hkInd'] = 1
  spkfData['zkInd'] = 2
  spkfData['xhat'] = np.array([ir0,hk0])                    # initial state
  spkfData['xhat'] = np.append(spkfData['xhat'],SOC0)
 
  # Covariance values
  spkfData['SigmaX'] = SigmaX0
  spkfData['SigmaV'] = SigmaV
  spkfData['SigmaW'] = SigmaW
  spkfData['Snoise'] = np.linalg.cholesky(np.diag((SigmaW,SigmaV))).real
  spkfData['Qbump'] = 5
  
#    SPKF specific parameters
  Nx = len(spkfData['xhat']) 
  spkfData['Nx'] = Nx                                        # state-vector length
  Ny = 1 
  Nu = 1 
  spkfData['Ny'] = Ny                                        # measurement-vector length
  spkfData['Nu'] = Nu                                        # input-vector length
  Nw = np.size(SigmaW)
  Nv = np.size(SigmaV)
  spkfData['Nw'] = Nw                                        # process-noise-vector length 
  spkfData['Nv'] = Nv                                        # sensor-noise-vector length
  Na = Nx+Nw+Nv 
  spkfData['Na'] = Na                                        # augmented-state-vector length
  
  h = np.sqrt(3)
  h = 3
  spkfData['h'] = h                                          # SPKF/CDKF tuning factor  
  Weight1 = (h*h-Na)/(h*h)                                   # weighting factors when computing mean and covariance
  Weight2 = (1/(2*h*h))
  set=np.ones((2*Na,1))
  spkfData['Wm'] = np.append(Weight1,Weight2*set)            # mean
  spkfData['Wc'] = spkfData['Wm']                            # covariance

  # previous value of current
  spkfData['priorI'] = 0
  spkfData['signIk'] = 0
  
  # store model data structure too
  spkfData['model'] = model
  return spkfData



