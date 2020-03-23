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
dataphugoidNMnotime = pd.read_csv('phugoidNM.csv', skiprows = 1, sep=',', names= ['dV_TAS','alpha','theta','q']) 
dataphugoidED = pd.read_csv('phugoidED.csv', skiprows = 1, sep=',', names= ['time','vane_AoA','Ahrs1_Pitch','Ahrs1_bPitchRate']) 

datashortperiodNMnotime = pd.read_csv('shortperiodNM.csv', skiprows = 1, sep=',', names= ['dV_TAS','alpha','theta','q'])
datashortperiodED = pd.read_csv('shortperiodED.csv', skiprows = 1, sep=',', names= ['time', 'vane_AoA','Ahrs1_Pitch','Ahrs1_bPitchRate'])

dataaperrollNMnotime = pd.read_csv('aperiodicrollNM.csv ', skiprows = 1, sep=',', names= ['beta','phi','p','r'])
dataaperrollED = pd.read_csv('aperiodicrollED.csv ', skiprows = 1, sep=',', names= ['time', 'Ahrs1_Roll','Ahrs1_bRollRate','Ahrs1_bYawRate'])

datadutchrollNMnotime = pd.read_csv('dutchrollNM.csv ', skiprows = 1, sep=',', names= ['beta','phi','p','r'])
datadutchrollED = pd.read_csv('dutchrollED.csv ', skiprows = 1, sep=',', names= ['time', 'Ahrs1_Roll','Ahrs1_bRollRate','Ahrs1_bYawRate'])

datadutchrollYDNMnotime = pd.read_csv('dutchrollYDNM.csv ', skiprows = 1, sep=',', names= ['beta','phi','p','r'])
datadutchrollYDED = pd.read_csv('dutchrollYDED.csv ', skiprows = 1, sep=',', names= ['time', 'Ahrs1_Roll','Ahrs1_bRollRate','Ahrs1_bYawRate'])

dataspiralNMnotime = pd.read_csv('spiralrollNM.csv ', skiprows = 1, sep=',', names= ['beta','phi','p','r'])
dataspiralED = pd.read_csv('spiralrollED.csv ', skiprows = 1, sep=',', names= ['time','Ahrs1_Roll','Ahrs1_bRollRate','Ahrs1_bYawRate'])

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
dataphugoidNM.plot(x = 'time', y = 'theta')

datadutchrollED.plot(x = 'time', y = 'Ahrs1_bRollRate')
datadutchrollNM.plot(x = 'time', y = 'p')

datadutchrollED.plot(x = 'time', y = 'Ahrs1_bYawRate')
datadutchrollNM.plot(x = 'time', y = 'r')

datadutchrollYDED.plot(x = 'time', y = 'Ahrs1_bRollRate')
datadutchrollYDNM.plot(x = 'time', y = 'p')

datadutchrollYDED.plot(x = 'time', y = 'Ahrs1_bYawRate')
datadutchrollYDNM.plot(x = 'time', y = 'r')

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











# dataphugoidnum['V_TAS']

#phugoidnum = pd.DataFrame(data=None, index=None, columns=None, dtype=None, copy=False)

#nrows = 6587, index_col = 0,

#phugoidnum = pd.DataFrame(data=phugoidNM.csv, index=None, columns=0, dtype=None, copy=False)

#timecolumn = []
#for i in range(:
#    timecolumn.append(i * 0.1)

#.astype(float)