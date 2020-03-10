# -*- coding: utf-8 -*-
"""
Created on Mon Mar  9 15:03:39 2020

@author: chang
"""
import numpy as np
from math import pi
from main import data
from parameters import m, e, CD0, CLa, A
import matplotlib.pyplot as plt

TAS = data.Dadc1_tas # true air speed
time = data.time
vane_AOA = data.vane_AOA # angle of attack
PressureAltitude = data.Dadc1_alt # [feet]
Pressure = ((850 - 1013) / 4813 * PressureAltitude + 1013) * 100 # Conversion from PA to P [Pa]

def CLCD_plot():
    CL  = np.linspace(0,1.7,500)
    CD = CD0 + CL**2 / (pi * A * e)
    
    plt.figure()
    plt.ylabel("CL")
    plt.xlabel("CD")
    plt.plot(CD,CL)
    plt.show()
    
CLCD_plot()
    
