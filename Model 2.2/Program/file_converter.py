from scipy.io import loadmat
import pandas as pd

matdata=loadmat(r'D:\Projects\Work\BMS\Models\Model 2\readonly\PANdata_P25.mat',simplify_cells=True)
modeldata=loadmat('PANmodel.mat', simplify_cells=True)
rawdata = matdata['DYNData']

eta = rawdata['eta']
Q=rawdata['Q']
time = rawdata['script1']['time']
deltat = time[1]-time[0]
time    = time-time[0]                                      # start time at 0
step = rawdata['script1']['step']
chgAh = rawdata['script1']['chgAh']
disAh = rawdata['script1']['disAh']   
current = rawdata['script1']['current']                     # discharge > 0; charge < 0.
voltage = rawdata['script1']['voltage']
rawTime = rawdata['script1']['rawTime']
rawCurrent = rawdata['script1']['rawCurrent']  
rawVoltage = rawdata['script1']['rawVoltage']
rawStep = rawdata['script1']['rawStep']
rawChgAh = rawdata['script1']['rawChgAh']
rawDisAh = rawdata['script1']['rawDisAh']    
soc     = rawdata['script1']['soc']

eta=pd.Series(eta)
Q=pd.Series(Q)
time=pd.Series(time)
deltat=pd.Series(deltat)
step=pd.Series(step)
chgAh=pd.Series(chgAh)
disAh=pd.Series(disAh)
current=pd.Series(current)
voltage=pd.Series(voltage)
rawTime=pd.Series(rawTime)
rawCurrent=pd.Series(rawCurrent)
rawVoltage=pd.Series(rawVoltage)
rawStep=pd.Series(rawStep)
rawChgAh=pd.Series(rawChgAh)
rawDisAh=pd.Series(rawDisAh)
soc=pd.Series(soc)

data=pd.concat([time,step,chgAh,disAh,current,voltage,rawTime,rawCurrent,rawVoltage,rawChgAh,rawDisAh,soc],
               keys=['time','step','chgAh','disAh','current','voltage','rawTime','rawCurrent','rawVoltage','rawChgAh','rawDisAh','soc'],axis=1)
script1=pd.DataFrame(data)
file=pd.concat([eta,Q],keys=['eta','Q'],axis=1)

with pd.ExcelWriter(r"D:\Projects\Work\BMS\Models\Model 2\readonly\PANdata_P25.xlsx") as writer:
        script1.to_excel(writer, sheet_name="script1", index=False)
        file.to_excel(writer, sheet_name="main", index=False)
        #data_frame3.to_excel(writer, sheet_name="Baked Items", index=False)
        
print("zhala")
    
###################################### model #######################################3

OCV0=modeldata['model']['OCV0']
OCVrel=modeldata['model']['OCVrel']
SOC=modeldata['model']['SOC']
OCV=modeldata['model']['OCV']
SOC0=modeldata['model']['SOC0']
SOCrel=modeldata['model']['SOCrel']
OCVeta=modeldata['model']['OCVeta']
OCVQ=modeldata['model']['OCVQ']
name=modeldata['model']['name']
temps=modeldata['model']['temps']
etaParam=modeldata['model']['etaParam']
QParam=modeldata['model']['QParam']
GParam=modeldata['model']['GParam']
M0Param=modeldata['model']['M0Param']
MParam=modeldata['model']['MParam']
R0Param=modeldata['model']['R0Param']
RCParam=modeldata['model']['RCParam']
RParam=modeldata['model']['RParam']
dOCV0=modeldata['model']['dOCV0']
dOCVrel=modeldata['model']['dOCVrel']

OCV0=pd.Series(OCV0)
OCVrel=pd.Series(OCVrel)
SOC=pd.Series(SOC)
OCV=pd.Series(OCV)
SOC0=pd.Series(SOC0)
SOCrel=pd.Series(SOCrel)
OCVeta=pd.Series(OCVeta)
OCVQ=pd.Series(OCVQ)
name=pd.Series(name)
temps=pd.Series(temps)
etaParam=pd.Series(etaParam)
QParam=pd.Series(QParam)
GParam=pd.Series(GParam)
M0Param=pd.Series(M0Param)
MParam=pd.Series(MParam)
R0Param=pd.Series(R0Param)
RCParam=pd.Series(RCParam)
RParam=pd.Series(RParam)
dOCV0=pd.Series(dOCV0)
dOCVrel=pd.Series(dOCVrel)

data2=pd.concat([OCV0,OCVrel,SOC,OCV,SOC0,SOCrel,OCVeta,OCVQ,name,temps,etaParam,QParam,GParam,M0Param,MParam,R0Param,RCParam,RParam,dOCV0,dOCVrel],
               keys=['OCV0','OCVrel','SOC','OCV','SOC0','SOCrel','OCVeta','OCVQ','name','temps','etaParam','QParam','GParam','M0Param','MParam','R0Param','RCParam','RParam','dOCV0','dOCVrel'],axis=1)
model=pd.DataFrame(data2)

with pd.ExcelWriter(r"D:\Projects\Work\BMS\Models\Model 2\readonly\PANmodel.xlsx") as writer:
         model.to_excel(writer, sheet_name="model", index=False)

print("Done")
# (sample syntax for writing to csv, excel file)
# file.to_csv(r'D:\Projects\Work\BMS\Models\Model 2\readonly\file.csv')
# file.to_excel(r"D:\Projects\Work\BMS\Models\Model 2\readonly\file.xlsx",sheet_name='script1')