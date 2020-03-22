# -*- coding: utf-8 -*-
"""
Created on Mon Mar  9 15:03:39 2020

@author: chang
"""
import numpy as np
import scipy as sp
import pandas as pd
import matplotlib.pyplot as plt
from main import m, e, CD0, AR, S, c
import matplotlib

# TAS = data.Dadc1_tas # true air speed
# time = data.time
# vane_AOA = data.vane_AOA # angle of attack
# PressureAltitude = data.Dadc1_alt # [feet]
# Pressure = ((850 - 1013) / 4813 * PressureAltitude + 1013) * 100 # Conversion from PA to P [Pa]



####### Data for Thrust.exe ########
# excel data
time_config1 = [16*60+32, 18*60+20, 20*60+6, 22*60, 24*60+48, 26*60+24] # [s]
IAS_config1 = np.array([248, 221, 188, 162, 140, 120]) * 0.514444 # [m/s]
pressure_alt = np.array([7000, 7000, 6980, 7000, 7000, 6980]) * 0.3048 # [m]
TAT_measured = np.array([1.7, 2.5, 3.7, 5.5, 7.7, 10.6])[::-1] + 273.15 # [Kelvin]

fuelflow_left = np.array([745, 641, 548, 456, 438, 456]) * 0.453592/60/60 # [kg/s]
fuelflow_right = np.array([803, 687, 593, 502, 472, 507]) * 0.453592/60/60 # [kg/s]

AOA = np.array([13.8, 11.8, 9.2, 7.8, 6.8, 5.8])[::-1]
AOA_rad = np.radians(AOA)

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

data_m = {'h_p': pressure_alt, 'Mach': Mach, 'temp diff': temp_diff, 'fuel flow left': fuelflow_left, 'fuel flow right': fuelflow_right }
matlab_data = pd.DataFrame(data_m)
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

Thrust = np.array(Thrust)        
V_TAS = Mach * 343
rho = 101325 / p0overp(pressure_alt) / 287 / TAT_corrected
Weight = (m - (fuelflow_left + fuelflow_right) * time_config1) * 9.80665

def Reynolds():
    # bulk_viscosity
    mu_0 = 1.7894e-5
    T_0 = 273.11 
    _S = 110.56
    mu = mu_0 *( TAT_corrected / T_0 ) ** (1.5) * (T_0 + _S) / (TAT_corrected + _S)

    return rho * V_TAS * c / mu 

Mach_min, Mach_max = np.min(Mach), np.max(Mach)
Reynolds_min, Reynolds_max = np.min(Reynolds()), np.max(Reynolds())

def CLCD_plot_stationary():
    
    CL = ( Weight - Thrust*np.sin(AOA_rad)) * 2 / (rho * V_TAS**2 * S)
    CD = ( Thrust * np.cos(AOA_rad)) * 2 / (rho * V_TAS**2 * S)
    
    ## Define text sizes for **SAVED** pictures (texpsize -- text export size)
    texpsize= [26,28,30]
    
    ## Input Arrays
    x = CD
    y = CL
    
    ## Best fit
    def func(CL, CD0, k):
        return CD0 + CL**2 * k
    
    cl_least = np.linspace(np.min(CL)-0.2, np.max(CL)+0.05,100)
    popt, pcov = sp.optimize.curve_fit(func, CL, CD)
    CD0 = popt[0]
    oswald = 1/ popt[1] / np.pi / AR
    
    # print(popt[0], popt[1])
    ## Graphing Parameters
    SMALL_SIZE  = texpsize[0]
    MEDIUM_SIZE = texpsize[1]
    BIGGER_SIZE = texpsize[2]
    
    plt.style.use('grayscale')
    plt.rc('font', size=MEDIUM_SIZE, family='serif')    ## controls default text sizes
    plt.rc('axes', titlesize=SMALL_SIZE)                ## fontsize of the axes title
    plt.rc('axes', labelsize=SMALL_SIZE)                ## fontsize of the x and y labels
    plt.rc('xtick', labelsize=SMALL_SIZE)               ## fontsize of the tick labels
    plt.rc('ytick', labelsize=SMALL_SIZE)               ## fontsize of the tick labels
    plt.rc('legend', fontsize=SMALL_SIZE)               ## legend fontsize
    plt.rc('figure', titlesize=BIGGER_SIZE)             ## fontsize of the figure title
    plt.rc('text', usetex=False)
    matplotlib.rcParams['lines.linewidth']  = 1.5
    matplotlib.rcParams['figure.facecolor'] = 'white'
    matplotlib.rcParams['axes.facecolor']   = 'white'
    matplotlib.rcParams["legend.fancybox"]  = False
    
    ## Graph
    fig, ax = plt.subplots(1,1,squeeze=False,figsize=(16,9))
    ax[0,0].scatter(x, y, label = 'Experimental data')
    ax[0,0].plot(func(cl_least, *popt), cl_least , label = 'Least Square regression $aCL^2 + b$')
    # ax[0,0].plot(label='Mach $\in (%1.3f, %1.3f)$ \n Reynolds $\in (%s, %s)$ '%(Mach_min, Mach_max, Reynolds_min, Reynolds_max))
    # ax[0,0].plot(x, y, label='(%f, %f)$ \n Reynolds (%f, %f)$ '%(Mach_min, Mach_max, Reynolds_min, Reynolds_max))
    # ax[0,0].set_title(r"Aircraft configuration: clean")
    # ax[0,0].plot(x+x, y+y, label="test2", linestyle="dashed")
    #ax[0,0].loglog(x, y, marker = "s", color='black', markerfacecolor='none', markeredgewidth=2, markersize=6, label="test")
    ax[0,0].set_ylabel(r"$C_L$")          ## String is treatable as latex code
    ax[0,0].set_xlabel(r"$C_D$")
    #ax[0,0].set_xlim(0,x[-1])
    ax[0,0].grid(True,which="major",color="#999999")
    ax[0,0].grid(True,which="minor",color="#DDDDDD",ls="--")
    ax[0,0].minorticks_on()
    ax[0,0].tick_params(which='major', length=10, width=2, direction='inout')
    ax[0,0].tick_params(which='minor', length=5, width=2, direction='in')
    ax[0,0].legend(loc=0, framealpha=1.0).get_frame().set_edgecolor('k')
    # ax[0,0].figtext(.8, .8, 'Mach $\in (%1.3f, %1.3f)$ \n Reynolds $\in (%s, %s)$' %(Mach_min, Mach_max, Reynolds_min, Reynolds_max))
    fig.savefig("CL-CD_cleanconfig.png", bbox_inches='tight')                                    ## Insert save destination
    
    ## If you want to see the figure, else disable last two lines.
    fig.tight_layout()
    plt.show()
    
    return CD0, oswald
    
def CLalpha_plot_stationary():
    
    CL = ( Weight - Thrust*np.sin(AOA_rad)) * 2 / (rho * V_TAS**2 * S)
    
    ## Define text sizes for **SAVED** pictures (texpsize -- text export size)
    texpsize= [26,28,30]
    
    ## Input Arrays
    x = AOA
    y = CL
    
    ## best fit
    def func(AOA, a, b):
        return a*AOA + b
        
    popt, pcov = sp.optimize.curve_fit(func, AOA, CL)
    CLalpha = popt[0]
    y_intercept = popt[1]
    aoa_least= np.linspace(np.min(AOA)-0.05, np.max(AOA)+0.1, 100)
    print(popt[0], popt[1])
    ## Graphing Parameters
    SMALL_SIZE  = texpsize[0]
    MEDIUM_SIZE = texpsize[1]
    BIGGER_SIZE = texpsize[2]
    
    plt.style.use('grayscale')
    plt.rc('font', size=MEDIUM_SIZE, family='serif')    ## controls default text sizes
    plt.rc('axes', titlesize=SMALL_SIZE)                ## fontsize of the axes title
    plt.rc('axes', labelsize=SMALL_SIZE)                ## fontsize of the x and y labels
    plt.rc('xtick', labelsize=SMALL_SIZE)               ## fontsize of the tick labels
    plt.rc('ytick', labelsize=SMALL_SIZE)               ## fontsize of the tick labels
    plt.rc('legend', fontsize=SMALL_SIZE)               ## legend fontsize
    plt.rc('figure', titlesize=BIGGER_SIZE)             ## fontsize of the figure title
    plt.rc('text', usetex=False)
    matplotlib.rcParams['lines.linewidth']  = 1.5
    matplotlib.rcParams['figure.facecolor'] = 'white'
    matplotlib.rcParams['axes.facecolor']   = 'white'
    matplotlib.rcParams["legend.fancybox"]  = False
    
    ## Graph
    fig, ax = plt.subplots(1,1,squeeze=False,figsize=(16,9))
    ax[0,0].scatter(x, y, label = 'Experimental data')
    ax[0,0].plot(aoa_least, func(aoa_least, *popt), label = 'Least Square regression $C_L = a\alpha + b$')
    # ax[0,0].set_title(r"Aircraft configuration: clean")
    # ax[0,0].plot(x,x*0.077-0.26, label="test2", linestyle="dashed")
    # ax[0,0].plot(x,x*0.07955175608200933-0.27201100205949574, label="test2", linestyle="dashed")
    #ax[0,0].loglog(x, y, marker = "s", color='black', markerfacecolor='none', markeredgewidth=2, markersize=6, label="test")
    ax[0,0].set_ylabel(r"$C_L$")          ## String is treatable as latex code
    ax[0,0].set_xlabel(r"$\alpha$ $\,\,[deg]$")
    #ax[0,0].set_xlim(0,x[-1])
    ax[0,0].grid(True,which="major",color="#999999")
    ax[0,0].grid(True,which="minor",color="#DDDDDD",ls="--")
    ax[0,0].minorticks_on()
    ax[0,0].tick_params(which='major', length=10, width=2, direction='inout')
    ax[0,0].tick_params(which='minor', length=5, width=2, direction='in')
    ax[0,0].legend(loc=0, framealpha=1.0).get_frame().set_edgecolor('k')
    fig.savefig("CL-alpha_cleanconfig.png", bbox_inches='tight')                                    ## Insert save destination
    
    ## If you want to see the figure, else disable last two lines.
    fig.tight_layout()
    plt.show()   
    
    return 0

CLCD_plot_stationary()
CLalpha_plot_stationary()


# Theoretical plot
def CL2CD_plot():
    CL  = np.linspace(0,1.0,100)
    CD = CLCD_plot_stationary()[0] + CL**2 / (np.pi * AR * CLCD_plot_stationary()[1])
    
    ## Define text sizes for **SAVED** pictures (texpsize -- text export size)
    texpsize= [26,28,30]
    
    ## Input Arrays
    x = np.linspace(1,10,10)
    y = np.ones((10))*2
    
    ## Graphing Parameters
    SMALL_SIZE  = texpsize[0]
    MEDIUM_SIZE = texpsize[1]
    BIGGER_SIZE = texpsize[2]
    
    plt.style.use('grayscale')
    plt.rc('font', size=MEDIUM_SIZE, family='serif')    ## controls default text sizes
    plt.rc('axes', titlesize=SMALL_SIZE)                ## fontsize of the axes title
    plt.rc('axes', labelsize=SMALL_SIZE)                ## fontsize of the x and y labels
    plt.rc('xtick', labelsize=SMALL_SIZE)               ## fontsize of the tick labels
    plt.rc('ytick', labelsize=SMALL_SIZE)               ## fontsize of the tick labels
    plt.rc('legend', fontsize=SMALL_SIZE)               ## legend fontsize
    plt.rc('figure', titlesize=BIGGER_SIZE)             ## fontsize of the figure title
    plt.rc('text', usetex=False)
    matplotlib.rcParams['lines.linewidth']  = 1.5
    matplotlib.rcParams['figure.facecolor'] = 'white'
    matplotlib.rcParams['axes.facecolor']   = 'white'
    matplotlib.rcParams["legend.fancybox"]  = False
    
    ## Graph
    fig, ax = plt.subplots(1,1,squeeze=False,figsize=(16,9))
    ax[0,0].plot(CD, CL)
    # ax[0,0].plot(x+x, y+y, label="test2", linestyle="dashed")
    #ax[0,0].loglog(x, y, marker = "s", color='black', markerfacecolor='none', markeredgewidth=2, markersize=6, label="test")
    ax[0,0].set_ylabel(r"$C_D$")          ## String is treatable as latex code
    ax[0,0].set_xlabel(r"$C_L$")
    #ax[0,0].set_xlim(0,x[-1])
    ax[0,0].grid(True,which="major",color="#999999")
    ax[0,0].grid(True,which="minor",color="#DDDDDD",ls="--")
    ax[0,0].minorticks_on()
    ax[0,0].tick_params(which='major', length=10, width=2, direction='inout')
    ax[0,0].tick_params(which='minor', length=5, width=2, direction='in')
    # ax[0,0].legend(loc=0, framealpha=1.0).get_frame().set_edgecolor('k')
    fig.savefig("CLCD-theoretical.png", bbox_inches='tight')                                    ## Insert save destination
    
    ## If you want to see the figure, else disable last two lines.
    fig.tight_layout()
    plt.show()            


CL2CD_plot()
    








