# -*- coding: utf-8 -*-
"""
Created on Mon Mar  9 15:03:39 2020

@author: chang
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# from main import data, m, e, CD0, CLa, A

# TAS = data.Dadc1_tas # true air speed
# time = data.time
# vane_AOA = data.vane_AOA # angle of attack
# PressureAltitude = data.Dadc1_alt # [feet]
# Pressure = ((850 - 1013) / 4813 * PressureAltitude + 1013) * 100 # Conversion from PA to P [Pa]

# def CLCD_plot():
#     CL  = np.linspace(0,1.7,500)
#     CD = CD0 + CL**2 / (np.pi * A * e)
    
#     plt.figure()
#     plt.ylabel("CL")
#     plt.xlabel("CD")
#     plt.plot(CD,CL)
#     plt.show()

####### Data for Thrust.exe ########
# excel data
time_config1 = [16*60+32, 18*60+20, 20*60+6, 22*60, 24*60+48, 26*60+24] # [s]
IAS_config1 = np.array([248, 221, 188, 162, 140, 120]) * 0.514444 # [m/s]
pressure_alt = np.array([7000, 7000, 6980, 7000, 7000, 6980]) * 0.3048 # [m]
TAT_measured = np.array([1.7, 2.5, 3.7, 5.5, 7.7, 10.6]) + 273.15 # [Kelvin]

fuelflow_left = np.array([745, 641, 548, 456, 438, 456]) * 0.453592/60/60 # [kg/s]
fuelflow_right = np.array([803, 687, 593, 502, 472, 507]) * 0.453592/60/60 # [kg/s]

_lambda = -0.0065
T_0 = 288.15
gamma = 1.4

def p0overp(h_p):
    g0 = 9.80665
    R = 287
    return 1 / (1 + _lambda * h_p / T_0)**(-1 * g0 /_lambda / R)

def Mach(p0p, IAS):
    subeqn = (1 + (gamma - 1) / 2 / gamma * 1.225 / 101325 * IAS ** 2) ** (gamma / (gamma - 1)) 
    Mach = np.sqrt( 2 / (gamma -1) * ( ((1 + p0p * (subeqn - 1)) ** ((gamma - 1) / gamma)) - 1) )
    return Mach

Mach = Mach(p0overp(pressure_alt), IAS_config1)

T_ISA = T_0 + _lambda * pressure_alt
TAT_corrected = TAT_measured / (1 + (gamma - 1) / 2 * Mach**2)
temp_diff = TAT_corrected - T_ISA
# temp_diff = T_ISA - TAT_corrected

data = {'h_p': pressure_alt, 'Mach': Mach, 'temp diff': temp_diff, 'fuel flow left': fuelflow_left, 'fuel flow right': fuelflow_right }
matlab_data = pd.DataFrame(data)
matlab_data.to_csv('thrust.exe/matlab.dat', header = False, index = False, sep = ' ')

Thrust = []
thrust_data = open("thrust.exe/thrust.dat")
lines = thrust_data.readlines()
for line in lines:
    separate_thrust = []
    for i in line.split():
        separate_thrust.append(float(i))
    thrust_sum = np.sum(separate_thrust)
    Thrust.append(thrust_sum)
        


# def CLCD_plot_stationary():
#     W = m * 9.80665
#     T = 
    
#     CL = ( W - T*np.sin(alpha)) * 2 / (rho, vel**2, S)
    
#     CD = ( T*np.cos(alpha)) * 2 / (rho, vel**2, S)
    
# CLCD_plot()
    
