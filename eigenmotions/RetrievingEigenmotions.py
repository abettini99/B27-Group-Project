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
dataphugoidNMnotime = pd.read_csv('flighttest_phugoidNM.csv', skiprows = 1, sep=',', names= ['dV_TAS','alpha','theta','q']) 
dataphugoidED = pd.read_csv('flighttest_phugoidED.csv', skiprows = 1, sep=',', names= ['time','vane_AoA','Ahrs1_Pitch','Ahrs1_bPitchRate']) 

datashortperiodNMnotime = pd.read_csv('flighttest_shortperiodNM.csv', skiprows = 1, sep=',', names= ['dV_TAS','alpha','theta','q'])
datashortperiodED = pd.read_csv('flighttest_shortperiodED.csv', skiprows = 1, sep=',', names= ['time', 'vane_AoA','Ahrs1_Pitch','Ahrs1_bPitchRate'])

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


dataphugoidED.plot(x = 'time', y = 'Ahrs1_Pitch')
#dataphugoidNM.plot(x = 'time', y = 'theta')
#
#datadutchrollED.plot(x = 'time', y = 'Ahrs1_bRollRate')
#datadutchrollNM.plot(x = 'time', y = 'p')
#
#datadutchrollED.plot(x = 'time', y = 'Ahrs1_bYawRate')
#datadutchrollNM.plot(x = 'time', y = 'r')
#
#datadutchrollYDED.plot(x = 'time', y = 'Ahrs1_bRollRate')
#datadutchrollYDNM.plot(x = 'time', y = 'p')
#
#datadutchrollYDED.plot(x = 'time', y = 'Ahrs1_bYawRate')
#datadutchrollYDNM.plot(x = 'time', y = 'r')

t = np.array(dataphugoidED['time'])


func = lambda t, lambda_Re1, lambda_Im1, lambda_Re2, lambda_Im2, a1, b1, a2, b2: np.exp(lambda_Re1 * t) * (a1*np.cos(lambda_Im1 * t) + b1*np.sin(lambda_Im1 * t)) + np.exp(lambda_Re2 * t) * (a2*np.cos(lambda_Im2 * t) + b2*np.sin(lambda_Im2 * t)) #+ np.exp(lambda_Re3 * t) * (a3*np.cos(lambda_Im3 * t) + b3*np.sin(lambda_Im3 * t)) + np.exp(lambda_Re4 * t) * (a4*np.cos(lambda_Im4 * t) + b4*np.sin(lambda_Im4 * t))

try:
    Output, CoVarMat = opt.curve_fit(func, dataphugoidED['time'], dataphugoidED['Ahrs1_Pitch'], bounds=[(-5, -5, -5, -5, -1, -1, -1, -1),(5, 5, 5, 5, 1, 1, 1, 1)])
except RuntimeError:
    print('ýo waddup homie')

#lambda_Re3, lambda_Im3, lambda_Re4, lambda_Im4,
#, a3, b3, a4, b4
#-5, -5, -5, -5, 
#, -1, -1, -1, -1
#5, 5, 5, 5, 1, 1, 1, 1,

#lambda_Re1 = Output[0]
#lambda_Im1 = Output[1]
#lambda_Re2 = Output[2]
#lambda_Im2 = Output[3]
#lambda_Re3 = Output[4]
#lambda_Im3 = Output[5]
#lambda_Re4 = Output[6]
#lambda_Im4 = Output[7]
#a1 = Output[8]
#b1 = Output[9]
#a2 = Output[10]
#b2 = Output[11]
#a3 = Output[12]
#b3 = Output[13]
#a4 = Output[14]
#b4 = Output[15]

lambda_Re1 = Output[0]
lambda_Im1 = Output[1]
lambda_Re2 = Output[2]
lambda_Im2 = Output[3]
a1 = Output[4]
b1 = Output[5]
a2 = Output[6]
b2 = Output[7]





#data = np.exp(lambda_Re1 * t) * (np.cos(lambda_Im1 * t) + np.sin(lambda_Im1 * t)) + np.exp(lambda_Re2 * t) * (np.cos(lambda_Im2 * t) + np.sin(lambda_Im2 * t)) + np.exp(lambda_Re3 * t) * (np.cos(lambda_Im3 * t) + np.sin(lambda_Im3 * t)) + np.exp(lambda_Re4 * t) * (np.cos(lambda_Im4 * t) + np.sin(lambda_Im4 * t))
data = ( np.exp(lambda_Re1 * t) * (a1*np.cos(lambda_Im1 * t) + b1*np.sin(lambda_Im1 * t)) 
   + np.exp(lambda_Re2 * t) * (a2*np.cos(lambda_Im2 * t) + b2*np.sin(lambda_Im2 * t)) )
    #+ np.exp(lambda_Re3 * t) * (a3*np.cos(lambda_Im3 * t) + b3*np.sin(lambda_Im3 * t)) 
    #+ np.exp(lambda_Re4 * t) * (a4*np.cos(lambda_Im4 * t) + b4*np.sin(lambda_Im4 * t)) )

#plt.figure()
plt.plot(t, data)
plt.show()




"""
phugoid: time vs theta
dutchroll: time vs p & time vs r
dutchrollYD: time vs p & time vs r
"""








#def func(CL, CD0, k):
#        return CD0 + CL**2 * k
#    
##    cl_least = np.linspace(np.min(CL)-0.23, np.max(CL)+0.1,100)
#    popt, pcov = sp.optimize.curve_fit(func, CL, CD)
#    CD0 = popt[0]
#    oswald = 1/ popt[1] / np.pi / AR


#error_func = lambda dt, K, beta : K * dt**beta
#ESVar, ESConv = opt.curve_fit(error_func, dataphugoidED['time'], dataphugoidED['vane_AoA'], bounds=[(0,0),(100,3)]) 



#func = lambda t, lambda_Re, lambda_Im, a, b : np.exp(lambda_Re * t) * (a*np.cos(lambda_Im * t) + b*np.sin(lambda_Im * t)) 
#Output, CoVarMat = opt.curve_fit(func, dataphugoidED['time'], dataphugoidED['vane_AoA'])

#EWVar, EWConv = opt.curve_fit(error_func, np.array([1,2,4,8,16])*dt, WeakConv[0,3,:], bounds=[(0,0.5),(100,3)])


# white noise and R^2









#https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.drop.html

# dataphugoidnum['V_TAS']

#phugoidnum = pd.DataFrame(data=None, index=None, columns=None, dtype=None, copy=False)

#nrows = 6587, index_col = 0,

#phugoidnum = pd.DataFrame(data=phugoidNM.csv, index=None, columns=0, dtype=None, copy=False)

#timecolumn = []
#for i in range(:
#    timecolumn.append(i * 0.1)

#.astype(float)

#test = np.array(dataphugoidED['Ahrs1_Pitch'])
#test2 = np.array(datadutchrollED['Ahrs1_bRollRate'])

#startphugoidED = dataphugoidED.where('Ahrs1_Pitch' > 0.)
#dataphugoidED.drop(labels=list(np.), axis=0)