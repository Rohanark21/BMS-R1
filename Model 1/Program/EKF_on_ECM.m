addpath readonly
load readonly/E2model

% retreive the dynamic hysteresis "gamma" value at 25 degC
hystGamma = getParamESC('GParam',-15,model)

% retreive the OCV at 50% SOC and 25 degrees celsius
ocv = OCVfromSOCtemp(0,-15,model)

RCGamma = getParamESC('RCParam',-15,model)

RGamma = getParamESC('RParam',-15,model)

QGamma = getParamESC('QParam',-15,model)

R0Gamma = getParamESC('R0Param',-15,model)

MGamma = getParamESC('MParam',-15,model)

M0Gamma = getParamESC('M0Param',-15,model)

etaGamma = getParamESC('etaParam',-15,model)

GGamma = getParamESC('GParam',-15,model)

% Find dOCV /dz at SOC = z from {SOC ,OCV } data
SOC=0:0.01:1;
function dOCVz = dOCVfromSOC (SOC ,OCV ,z)
dZ = SOC (2) - SOC (1) ; % Find spacing of SOC vector
dUdZ = diff ( OCV )/dZ; % Scaled forward finite difference
dOCV = ([dUdZ(1) dUdZ] + [dUdZ dUdZ(end)]) /2; % Avg of fwd / bkwd diffs
dOCVz = interp1 (SOC ,dOCV ,z); % Could make more efficient than this .
end
