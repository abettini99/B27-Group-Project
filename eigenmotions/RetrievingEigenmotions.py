# -*- Retrieving Eigenmotions (through P and T_(1/2) -*-
"""
Created on Mon Mar 23 15:03:05 2020

@author: Natascha
"""
import pandas as pd
import matplotlib.pyplot as plt
import scipy.optimize as opt
import numpy as np
from math import *

#Reading the csv files
dataphugoidNMnotime = pd.read_csv('flighttest_phugoidNM.csv', skiprows = 1, sep=',', names= ['dV_TAS','alpha','theta','q']) #V_{TAS} - V_0
dataphugoidED = pd.read_csv('flighttest_phugoidED.csv', skiprows = 1, sep=',', names= ['time','Dadc1_tas', 'vane_AoA','Ahrs1_Pitch','Ahrs1_bPitchRate']) #time,Dadc1_tas,vane_AoA,Ahrs1_Pitch,Ahrs1_bPitchRate

datashortperiodNMnotime = pd.read_csv('flighttest_shortperiodNM.csv', skiprows = 1, sep=',', names= ['dV_TAS','alpha','theta','q'])
datashortperiodED = pd.read_csv('flighttest_shortperiodED.csv', skiprows = 1, sep=',', names= ['time', 'Dadc1_tas', 'vane_AoA','Ahrs1_Pitch','Ahrs1_bPitchRate']) #time,Dadc1_tas,vane_AoA,Ahrs1_Pitch,Ahrs1_bPitchRate

dataaperrollNMnotime = pd.read_csv('flighttest_aperiodicrollNM.csv ', skiprows = 1, sep=',', names= ['beta','phi','p','r'])
dataaperrollED = pd.read_csv('flighttest_aperiodicrollED.csv ', skiprows = 1, sep=',', names= ['time', 'Ahrs1_Roll','Ahrs1_bRollRate','Ahrs1_bYawRate']) 

datadutchrollNMnotime = pd.read_csv('flighttest_dutchrollNM.csv ', skiprows = 1, sep=',', names= ['beta','phi','p','r'])
datadutchrollED = pd.read_csv('flighttest_dutchrollED.csv ', skiprows = 1, sep=',', names= ['time', 'Ahrs1_Roll','Ahrs1_bRollRate','Ahrs1_bYawRate'])

datadutchrollYDNMnotime = pd.read_csv('flighttest_dutchrollYDNM.csv ', skiprows = 1, sep=',', names= ['beta','phi','p','r'])
datadutchrollYDED = pd.read_csv('flighttest_dutchrollYDED.csv ', skiprows = 1, sep=',', names= ['time', 'Ahrs1_Roll','Ahrs1_bRollRate','Ahrs1_bYawRate'])

dataspiralNMnotime = pd.read_csv('flighttest_spiralNM.csv ', skiprows = 1, sep=',', names= ['beta','phi','p','r'])
dataspiralED = pd.read_csv('flighttest_spiralED.csv ', skiprows = 1, sep=',', names= ['time','Ahrs1_Roll','Ahrs1_bRollRate','Ahrs1_bYawRate'])

#Adding time column to numerical data
phugoidtime = dataphugoidED['time']
shortperiodtime = datashortperiodED['time']
aperrolltime = dataaperrollED['time']
dutchrolltime = datadutchrollED['time']
dutchrollYDtime = datadutchrollYDED['time']
spiraltime = dataspiralED['time']

dataphugoidNM = dataphugoidNMnotime.join(phugoidtime)
datashortperiodNM = datashortperiodNMnotime.join(shortperiodtime)
dataaperrollNM = dataaperrollNMnotime.join(aperrolltime)
datadutchrollNM = datadutchrollNMnotime.join(dutchrolltime)
datadutchrollYDNM = datadutchrollYDNMnotime.join(dutchrollYDtime)
dataspiralNM = dataspiralNMnotime.join(spiraltime)

#Visualising the data
#dataphugoidNM.plot(x = 'time', y = 'dV_TAS')
#dataphugoidNM.plot(x = 'time', y = 'alpha')
#dataphugoidNM.plot(x = 'time', y = 'theta')
#dataphugoidNM.plot(x = 'time', y = 'q')
#
#dataphugoidED.plot(x = 'time', y = 'Dadc1_tas')
#dataphugoidED.plot(x = 'time', y = 'vane_AoA')
#dataphugoidED.plot(x = 'time', y = 'Ahrs1_Pitch')
#dataphugoidED.plot(x = 'time', y = 'Ahrs1_bPitchRate')
#
#
datashortperiodNM.plot(x = 'time', y = 'dV_TAS')
datashortperiodNM.plot(x = 'time', y = 'alpha')
datashortperiodNM.plot(x = 'time', y = 'theta')
datashortperiodNM.plot(x = 'time', y = 'q')
#
#datashortperiodED.plot(x = 'time', y = 'Dadc1_tas')
#datashortperiodED.plot(x = 'time', y = 'vane_AoA')
#datashortperiodED.plot(x = 'time', y = 'Ahrs1_Pitch')
#datashortperiodED.plot(x = 'time', y = 'Ahrs1_bPitchRate')
#
#
#dataaperrollNM.plot(x = 'time', y = 'beta')
#dataaperrollNM.plot(x = 'time', y = 'phi')
#dataaperrollNM.plot(x = 'time', y = 'p')
#dataaperrollNM.plot(x = 'time', y = 'r')
#
#dataaperrollED.plot(x = 'time', y = 'Ahrs1_Roll')
#dataaperrollED.plot(x = 'time', y = 'Ahrs1_bRollRate')
#dataaperrollED.plot(x = 'time', y = 'Ahrs1_bYawRate')
#
#
#datadutchrollNM.plot(x = 'time', y = 'beta')
#datadutchrollNM.plot(x = 'time', y = 'phi')
#datadutchrollNM.plot(x = 'time', y = 'p')
#datadutchrollNM.plot(x = 'time', y = 'r')
#
#datadutchrollED.plot(x = 'time', y = 'Ahrs1_Roll')
#datadutchrollED.plot(x = 'time', y = 'Ahrs1_bRollRate')
#datadutchrollED.plot(x = 'time', y = 'Ahrs1_bYawRate')
#
#
#datadutchrollYDNM.plot(x = 'time', y = 'beta')
#datadutchrollYDNM.plot(x = 'time', y = 'phi')
#datadutchrollYDNM.plot(x = 'time', y = 'p')
#datadutchrollYDNM.plot(x = 'time', y = 'r')
#
#datadutchrollYDED.plot(x = 'time', y = 'Ahrs1_Roll')
#datadutchrollYDED.plot(x = 'time', y = 'Ahrs1_bRollRate')
#datadutchrollYDED.plot(x = 'time', y = 'Ahrs1_bYawRate')
#
#
#dataspiralNM.plot(x = 'time', y = 'beta')
#dataspiralNM.plot(x = 'time', y = 'phi')
#dataspiralNM.plot(x = 'time', y = 'p')
#dataspiralNM.plot(x = 'time', y = 'r')
#
#dataspiralED.plot(x = 'time', y = 'Ahrs1_Roll')
#dataspiralED.plot(x = 'time', y = 'Ahrs1_bRollRate')
#dataspiralED.plot(x = 'time', y = 'Ahrs1_bYawRate')




#PhugoidNM
phugoidNMtime = list(dataphugoidNM['time'])
phugoidNMtheta = list(dataphugoidNM['theta'])
a0_phugoidNM = min(phugoidNMtheta)
t0_phugoidNM = phugoidNMtime[phugoidNMtheta.index(a0_phugoidNM)]

beginsliceval_phugoidNM = 2720
endsliceval_phugoidNM = 2760

index_beginslice_phugoidNM = phugoidNMtime.index(beginsliceval_phugoidNM)
index_endslice_phugoidNM = phugoidNMtime.index(endsliceval_phugoidNM)
slicedphugoidNMtheta = phugoidNMtheta[index_beginslice_phugoidNM:index_endslice_phugoidNM]
y_phugoidNM = min(slicedphugoidNMtheta)
t1_phugoidNM = phugoidNMtime[phugoidNMtheta.index(y_phugoidNM)]
t_phugoidNM = t1_phugoidNM - t0_phugoidNM

lambda_Re_phugoidNM = np.log(y_phugoidNM/a0_phugoidNM) / t_phugoidNM
lambda_Im_phugoidNM = 2 * np.pi / t_phugoidNM

#PhugoidED
phugoidEDtime = list(dataphugoidED['time'])
phugoidEDtheta = list(dataphugoidED['Ahrs1_Pitch'])
a0_phugoidED = min(phugoidEDtheta)
t0_phugoidED = phugoidEDtime[phugoidEDtheta.index(a0_phugoidED)]

beginsliceval_phugoidED = 2720
endsliceval_phugoidED = 2760

index_beginslice_phugoidED = phugoidEDtime.index(beginsliceval_phugoidED)
index_endslice_phugoidED = phugoidEDtime.index(endsliceval_phugoidED)
slicedphugoidEDtheta = phugoidEDtheta[index_beginslice_phugoidED:index_endslice_phugoidED]
y_phugoidED = min(slicedphugoidEDtheta)
t1_phugoidED = phugoidEDtime[phugoidEDtheta.index(y_phugoidED)]
t_phugoidED = t1_phugoidED - t0_phugoidED

lambda_Re_phugoidED = np.log(y_phugoidED/a0_phugoidED) / t_phugoidED
lambda_Im_phugoidED = 2 * np.pi / t_phugoidED

#Dutch roll NM
dutchrollNMtime = list(datadutchrollNM['time'])
dutchrollNMr = list(datadutchrollNM['r'])

beginsliceval1_dutchrollNM = 3026
endsliceval1_dutchrollNM = 3028

index_beginslice1_dutchrollNM = dutchrollNMtime.index(beginsliceval1_dutchrollNM)
index_endslice1_dutchrollNM = dutchrollNMtime.index(endsliceval1_dutchrollNM)
sliced1dutchrollNMr = dutchrollNMr[index_beginslice1_dutchrollNM:index_endslice1_dutchrollNM]

a0_dutchrollNM = max(sliced1dutchrollNMr)
t0_dutchrollNM = dutchrollNMtime[dutchrollNMr.index(a0_dutchrollNM)]

beginsliceval2_dutchrollNM = 3029
endsliceval2_dutchrollNM = 3031

index_beginslice2_dutchrollNM = dutchrollNMtime.index(beginsliceval2_dutchrollNM)
index_endslice2_dutchrollNM = dutchrollNMtime.index(endsliceval2_dutchrollNM)
sliced2dutchrollNMr = dutchrollNMr[index_beginslice2_dutchrollNM:index_endslice2_dutchrollNM]
y_dutchrollNM = max(sliced2dutchrollNMr)
t1_dutchrollNM = dutchrollNMtime[dutchrollNMr.index(y_dutchrollNM)]
t_dutchrollNM = t1_dutchrollNM - t0_dutchrollNM

lambda_Re_dutchrollNM = np.log(y_dutchrollNM/a0_dutchrollNM) / t_dutchrollNM
lambda_Im_dutchrollNM = 2 * np.pi / t_dutchrollNM

#Dutch roll ED
dutchrollEDtime = list(datadutchrollED['time'])
dutchrollEDr = list(datadutchrollED['Ahrs1_bYawRate'])

beginsliceval1_dutchrollED = 3026
endsliceval1_dutchrollED = 3028

index_beginslice1_dutchrollED = dutchrollEDtime.index(beginsliceval1_dutchrollED)
index_endslice1_dutchrollED = dutchrollEDtime.index(endsliceval1_dutchrollED)
sliced1dutchrollEDr = dutchrollEDr[index_beginslice1_dutchrollED:index_endslice1_dutchrollED]

a0_dutchrollED = max(sliced1dutchrollEDr)
t0_dutchrollED = dutchrollEDtime[dutchrollEDr.index(a0_dutchrollED)]

beginsliceval2_dutchrollED = 3030
endsliceval2_dutchrollED = 3031

index_beginslice2_dutchrollED = dutchrollEDtime.index(beginsliceval2_dutchrollED)
index_endslice2_dutchrollED = dutchrollEDtime.index(endsliceval2_dutchrollED)
sliced2dutchrollEDr = dutchrollEDr[index_beginslice2_dutchrollED:index_endslice2_dutchrollED]
y_dutchrollED = max(sliced2dutchrollEDr)
t1_dutchrollED = dutchrollEDtime[dutchrollEDr.index(y_dutchrollED)]
t_dutchrollED = t1_dutchrollED - t0_dutchrollED

lambda_Re_dutchrollED = np.log(y_dutchrollED/a0_dutchrollED) / t_dutchrollED
lambda_Im_dutchrollED = 2 * np.pi / t_dutchrollED

#Dutch roll YD NM
dutchrollYDNMtime = list(datadutchrollYDNM['time'])
dutchrollYDNMr = list(datadutchrollYDNM['r'])

beginsliceval1_dutchrollYDNM = 3095
endsliceval1_dutchrollYDNM = 3096

index_beginslice1_dutchrollYDNM = dutchrollYDNMtime.index(beginsliceval1_dutchrollYDNM)
index_endslice1_dutchrollYDNM = dutchrollYDNMtime.index(endsliceval1_dutchrollYDNM)
sliced1dutchrollYDNMr = dutchrollYDNMr[index_beginslice1_dutchrollYDNM:index_endslice1_dutchrollYDNM]

a0_dutchrollYDNM = max(sliced1dutchrollYDNMr)
t0_dutchrollYDNM = dutchrollYDNMtime[dutchrollYDNMr.index(a0_dutchrollYDNM)]

beginsliceval2_dutchrollYDNM = 3098
endsliceval2_dutchrollYDNM = 3099

index_beginslice2_dutchrollYDNM = dutchrollYDNMtime.index(beginsliceval2_dutchrollYDNM)
index_endslice2_dutchrollYDNM = dutchrollYDNMtime.index(endsliceval2_dutchrollYDNM)
sliced2dutchrollYDNMr = dutchrollYDNMr[index_beginslice2_dutchrollYDNM:index_endslice2_dutchrollYDNM]
y_dutchrollYDNM = max(sliced2dutchrollYDNMr)
t1_dutchrollYDNM = dutchrollYDNMtime[dutchrollYDNMr.index(y_dutchrollYDNM)]
t_dutchrollYDNM = t1_dutchrollYDNM - t0_dutchrollYDNM

lambda_Re_dutchrollYDNM = np.log(y_dutchrollYDNM/a0_dutchrollYDNM) / t_dutchrollYDNM
lambda_Im_dutchrollYDNM = 2 * np.pi / t_dutchrollYDNM

#Dutch roll YD ED
dutchrollYDEDtime = list(datadutchrollYDED['time'])
dutchrollYDEDr = list(datadutchrollYDED['Ahrs1_bYawRate'])

beginsliceval1_dutchrollYDED = 3095
endsliceval1_dutchrollYDED = 3098

index_beginslice1_dutchrollYDED = dutchrollYDEDtime.index(beginsliceval1_dutchrollYDED)
index_endslice1_dutchrollYDED = dutchrollYDEDtime.index(endsliceval1_dutchrollYDED)
sliced1dutchrollYDEDr = dutchrollYDEDr[index_beginslice1_dutchrollYDED:index_endslice1_dutchrollYDED]

a0_dutchrollYDED = max(sliced1dutchrollYDEDr)
t0_dutchrollYDED = dutchrollYDEDtime[dutchrollYDEDr.index(a0_dutchrollYDED)]

beginsliceval2_dutchrollYDED = 3098
endsliceval2_dutchrollYDED = 3100

index_beginslice2_dutchrollYDED = dutchrollYDEDtime.index(beginsliceval2_dutchrollYDED)
index_endslice2_dutchrollYDED = dutchrollYDEDtime.index(endsliceval2_dutchrollYDED)
sliced2dutchrollYDEDr = dutchrollYDEDr[index_beginslice2_dutchrollYDED:index_endslice2_dutchrollYDED]
y_dutchrollYDED = max(sliced2dutchrollYDEDr)
t1_dutchrollYDED = dutchrollYDEDtime[dutchrollYDEDr.index(y_dutchrollYDED)]
t_dutchrollYDED = t1_dutchrollYDED - t0_dutchrollYDED

lambda_Re_dutchrollYDED = np.log(y_dutchrollYDED/a0_dutchrollYDED) / t_dutchrollYDED
lambda_Im_dutchrollYDED = 2 * np.pi / t_dutchrollYDED

#Short period NM
shortperiodNMtime = list(datashortperiodNM['time'])
shortperiodNMq = list(datashortperiodNM['q'])

beginsettime_shortperiodNM = 2632
endsettime_shortperiodNM = 2634

index_beginsettime_shortperiodNM = shortperiodNMtime.index(beginsettime_shortperiodNM)
index_endsettime_shortperiodNM = shortperiodNMtime.index(endsettime_shortperiodNM)

a_shortperiodNM = shortperiodNMq[index_beginsettime_shortperiodNM]
y_shortperiodNM = shortperiodNMq[index_endsettime_shortperiodNM]

lambda_shortperiodNM = np.log(y_shortperiodNM / a_shortperiodNM) / (endsettime_shortperiodNM - beginsettime_shortperiodNM)

#Short period ED
shortperiodEDtime = list(datashortperiodED['time'])
shortperiodEDq = list(datashortperiodED['Ahrs1_bPitchRate'])

beginsettime_shortperiodED = 2632
endsettime_shortperiodED = 2634

index_beginsettime_shortperiodED = shortperiodEDtime.index(beginsettime_shortperiodED)
index_endsettime_shortperiodED = shortperiodEDtime.index(endsettime_shortperiodED)

a_shortperiodED = shortperiodEDq[index_beginsettime_shortperiodED]
y_shortperiodED = shortperiodEDq[index_endsettime_shortperiodED]

lambda_shortperiodED = np.log(y_shortperiodED / a_shortperiodED) / (endsettime_shortperiodED - beginsettime_shortperiodED)

#Aperiodic roll NM
aperrollNMtime = list(dataaperrollNM['time'])
aperrollNMphi = list(dataaperrollNM['phi'])

beginsettime_aperrollNM = 2904
endsettime_aperrollNM = 2908

index_beginsettime_aperrollNM = aperrollNMtime.index(beginsettime_aperrollNM)
index_endsettime_aperrollNM = aperrollNMtime.index(endsettime_aperrollNM)

a_aperrollNM = aperrollNMphi[index_beginsettime_aperrollNM]
y_aperrollNM = aperrollNMphi[index_endsettime_aperrollNM]

lambda_aperrollNM = np.log(y_aperrollNM / a_aperrollNM) / (endsettime_aperrollNM - beginsettime_aperrollNM)

#Aperiodic roll ED
aperrollEDtime = list(dataaperrollED['time'])
aperrollEDphi = list(dataaperrollED['Ahrs1_Roll'])

beginsettime_aperrollED = 2905.5
endsettime_aperrollED = 2908

index_beginsettime_aperrollED = aperrollEDtime.index(beginsettime_aperrollED)
index_endsettime_aperrollED = aperrollEDtime.index(endsettime_aperrollED)

a_aperrollED = aperrollEDphi[index_beginsettime_aperrollED]
y_aperrollED = aperrollEDphi[index_endsettime_aperrollED]

lambda_aperrollED = np.log(y_aperrollED / a_aperrollED) / (endsettime_aperrollED - beginsettime_aperrollED)

#Spiral NM 
spiralNMtime = list(dataspiralNM['time'])
spiralNMphi = list(dataspiralNM['phi'])

beginsettime_spiralNM = 3310
endsettime_spiralNM = 3380

index_beginsettime_spiralNM = spiralNMtime.index(beginsettime_spiralNM)
index_endsettime_spiralNM = spiralNMtime.index(endsettime_spiralNM)

a_spiralNM = spiralNMphi[index_beginsettime_spiralNM]
y_spiralNM = spiralNMphi[index_endsettime_spiralNM]

lambda_spiralNM = np.log(y_spiralNM / a_spiralNM) / (endsettime_spiralNM - beginsettime_spiralNM)

#Spiral ED
spiralEDtime = list(dataspiralED['time'])
spiralEDphi = list(dataspiralED['Ahrs1_Roll'])

beginsettime_spiralED = 3310
endsettime_spiralED = 3380

index_beginsettime_spiralED = spiralEDtime.index(beginsettime_spiralED)
index_endsettime_spiralED = spiralEDtime.index(endsettime_spiralED)

a_spiralED = spiralEDphi[index_beginsettime_spiralED]
y_spiralED = spiralEDphi[index_endsettime_spiralED]

lambda_spiralED = np.log(y_spiralED / a_spiralED) / (endsettime_spiralED - beginsettime_spiralED)

# =============================================================================
# t = np.array(dataphugoidED['time'])[::50]
# y = dataphugoidED['Ahrs1_Pitch']
# 
# plt.plot
# 
# func = lambda t, lambda_Re1, lambda_Im1, lambda_Re2, lambda_Im2, a1, b1, a2, b2: np.exp(lambda_Re1 * t) * (a1*np.cos(lambda_Im1 * t) + b1*np.sin(lambda_Im1 * t)) + np.exp(lambda_Re2 * t) * (a2*np.cos(lambda_Im2 * t) + b2*np.sin(lambda_Im2 * t)) #+ np.exp(lambda_Re3 * t) * (a3*np.cos(lambda_Im3 * t) + b3*np.sin(lambda_Im3 * t)) + np.exp(lambda_Re4 * t) * (a4*np.cos(lambda_Im4 * t) + b4*np.sin(lambda_Im4 * t))
# 
# try:
#     Output, CoVarMat = opt.curve_fit(func, dataphugoidED['time'], dataphugoidED['Ahrs1_Pitch'], bounds=[(-5, -5, -5, -5, -1, -1, -1, -1),(5, 5, 5, 5, 1, 1, 1, 1)])
# except RuntimeError:
#     print('ýo waddup homie')
# 
# #lambda_Re3, lambda_Im3, lambda_Re4, lambda_Im4,
# #, a3, b3, a4, b4
# #-5, -5, -5, -5, 
# #, -1, -1, -1, -1
# #5, 5, 5, 5, 1, 1, 1, 1,
# 
# #lambda_Re1 = Output[0]
# #lambda_Im1 = Output[1]
# #lambda_Re2 = Output[2]
# #lambda_Im2 = Output[3]
# #lambda_Re3 = Output[4]
# #lambda_Im3 = Output[5]
# #lambda_Re4 = Output[6]
# #lambda_Im4 = Output[7]
# #a1 = Output[8]
# #b1 = Output[9]
# #a2 = Output[10]
# #b2 = Output[11]
# #a3 = Output[12]
# #b3 = Output[13]
# #a4 = Output[14]
# #b4 = Output[15]
# 
# lambda_Re1 = Output[0]
# lambda_Im1 = Output[1]
# lambda_Re2 = Output[2]
# lambda_Im2 = Output[3]
# a1 = Output[4]
# b1 = Output[5]
# a2 = Output[6]
# b2 = Output[7]
# 
# 
# 
# 
# 
# #data = np.exp(lambda_Re1 * t) * (np.cos(lambda_Im1 * t) + np.sin(lambda_Im1 * t)) + np.exp(lambda_Re2 * t) * (np.cos(lambda_Im2 * t) + np.sin(lambda_Im2 * t)) + np.exp(lambda_Re3 * t) * (np.cos(lambda_Im3 * t) + np.sin(lambda_Im3 * t)) + np.exp(lambda_Re4 * t) * (np.cos(lambda_Im4 * t) + np.sin(lambda_Im4 * t))
# data = ( np.exp(lambda_Re1 * t) * (a1*np.cos(lambda_Im1 * t) + b1*np.sin(lambda_Im1 * t)) 
#    + np.exp(lambda_Re2 * t) * (a2*np.cos(lambda_Im2 * t) + b2*np.sin(lambda_Im2 * t)) )
#     #+ np.exp(lambda_Re3 * t) * (a3*np.cos(lambda_Im3 * t) + b3*np.sin(lambda_Im3 * t)) 
#     #+ np.exp(lambda_Re4 * t) * (a4*np.cos(lambda_Im4 * t) + b4*np.sin(lambda_Im4 * t)) )
# 
# #plt.figure()
# plt.plot(t, data)
# plt.show()
# =============================================================================



# =============================================================================
# #func = lambda t, lambda_Re, lambda_Im, a, b : np.exp(lambda_Re * t) * (a*np.cos(lambda_Im * t) + b*np.sin(lambda_Im * t)) 
# #Output, CoVarMat = opt.curve_fit(func, dataphugoidED['time'], dataphugoidED['vane_AoA'])
# 
# #EWVar, EWConv = opt.curve_fit(error_func, np.array([1,2,4,8,16])*dt, WeakConv[0,3,:], bounds=[(0,0.5),(100,3)])
# 
# #dataphugoidNM.loc(dataphugoidNM['theta'] == min(dataphugoidNM['theta']))
#
# 
# # white noise and R^2
# 
# 
# #https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.drop.html
# 
# 
# #timecolumn = []
# #for i in range(:
# #    timecolumn.append(i * 0.1)
# 
# #.astype(float)
# 
# #test = np.array(dataphugoidED['Ahrs1_Pitch'])
# #test2 = np.array(datadutchrollED['Ahrs1_bRollRate'])
# 
# #startphugoidED = dataphugoidED.where('Ahrs1_Pitch' > 0.)
# #dataphugoidED.drop(labels=list(np.), axis=0)
# =============================================================================

#index_roots_phugoidNM = []
#roots_phugoidNM = []
#roots_time_phugoidNM = []
#for i in phugoidNMtheta:
#    if np.abs(i) < 8.0e-4:
#        index_roots_phugoidNM.append(phugoidNMtheta.index(i))
#        roots_phugoidNM.append(i)
#        roots_time_phugoidNM.append(phugoidNMtime[phugoidNMtheta.index(i)])
##        print(i)
##plt.figure()
##plt.plot(phugoidNMtime, phugoidNMtheta)
##plt.scatter(roots_time_phugoid, roots_phugoid)
##plt.show()
#        
#P_phugoidNM = roots_time_phugoidNM[9]-roots_time_phugoidNM[7]