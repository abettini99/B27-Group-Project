# -*- Retrieving Eigenmotions (through P and T_(1/2) -*-
"""
Created on Mon Mar 23 15:03:05 2020

@author: Natascha
"""
import pandas as pd
import matplotlib.pyplot as plt

#timecolumn = []
#for i in range(:
#    timecolumn.append(i * 0.1)


dataphugoidNM = pd.read_csv('phugoidNM.csv', skiprows = 0, sep=',', names= ['dV_TAS','alpha','theta','q'])
dataphugoidED = pd.read_csv('phugoidED.csv', skiprows = 0, sep=',', names= ['vane_AoA','Ahrs1_Pitch','Ahrs1_bPitchRate'])

datashortperiodNM = pd.read_csv('shortperiodNM.csv', skiprows = 0, sep=',', names= ['dV_TAS','alpha','theta','q'])
datashortperiodED = pd.read_csv('shortperiodED.csv', skiprows = 0, sep=',', names= ['vane_AoA','Ahrs1_Pitch','Ahrs1_bPitchRate'])

dataaperrollNM = pd.read_csv('aperiodicrollNM.csv ', skiprows = 0, sep=',', names= ['beta','phi','p','r'])
dataaperrollED = pd.read_csv('aperiodicrollED.csv ', skiprows = 0, sep=',', names= ['Ahrs1_Roll','Ahrs1_bRollRate','Ahrs1_bYawRate'])

datadutchrollNM = pd.read_csv('dutchrollNM.csv ', skiprows = 0, sep=',', names= ['beta','phi','p','r'])
datadutchrollED = pd.read_csv('dutchrollED.csv ', skiprows = 0, sep=',', names= ['Ahrs1_Roll','Ahrs1_bRollRate','Ahrs1_bYawRate'])

datadutchrollYDNM = pd.read_csv('dutchrollYDNM.csv ', skiprows = 0, sep=',', names= ['beta','phi','p','r'])
datadutchrollYDED = pd.read_csv('dutchrollYDED.csv ', skiprows = 0, sep=',', names= ['Ahrs1_Roll','Ahrs1_bRollRate','Ahrs1_bYawRate'])

dataspiralNM = pd.read_csv('spiralrollYDNM.csv ', skiprows = 0, sep=',', names= ['beta','phi','p','r'])
dataspiralED = pd.read_csv('spiralrollYDED.csv ', skiprows = 0, sep=',', names= ['Ahrs1_Roll','Ahrs1_bRollRate','Ahrs1_bYawRate'])

#dataphugoidNM_new = manouvre(data, 'aperroll')
dataphugoidNM.plot(x='time', y='Ahrs1_bRollRate')














# dataphugoidnum['V_TAS']

#phugoidnum = pd.DataFrame(data=None, index=None, columns=None, dtype=None, copy=False)

#nrows = 6587, index_col = 0,

#phugoidnum = pd.DataFrame(data=phugoidNM.csv, index=None, columns=0, dtype=None, copy=False)