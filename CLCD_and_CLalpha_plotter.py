"""
Created on Mon Mar  9 15:03:39 2020

@author: chang
"""
import numpy as np
import scipy as sp
import scipy.optimize as spopt
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib

# TAS = data.Dadc1_tas # true air speed
# time = data.time
# vane_AOA = data.vane_AOA # angle of attack
# PressureAltitude = data.Dadc1_alt # [feet]
# Pressure = ((850 - 1013) / 4813 * PressureAltitude + 1013) * 100 # Conversion from PA to P [Pa]
b = 15.911
S = 30.00
AR = b**2 / S
c = 2.0569

####### Data for Thrust.exe ########
# experimental data
print("Experiemental")
m = 9165 * 0.453592 + 90 + 102 + 80 + 83 + 94 + 84 + 74 + 79 + 103 + 4100 * 0.453592 # [kg]
time_config1 = [16*60+32, 18*60+20, 20*60+6, 22*60, 24*60+48, 26*60+24] # [s]
IAS_config1 = (np.array([248, 221, 188, 162, 140, 120]) - 2) * 0.514444 # [m/s]
pressure_alt = np.array([7000, 7000, 6980, 7000, 7000, 6980]) * 0.3048 # [m]
AOA = np.array([1.7, 2.5, 3.7, 5.5, 7.7, 10.6])

fuelflow_left = np.array([745, 641, 548, 456, 438, 456]) * 0.453592/60/60 # [kg/s]
fuelflow_right = np.array([803, 687, 593, 502, 472, 507]) * 0.453592/60/60 # [kg/s]
fuelused = np.array([367, 400, 430, 470, 497, 515]) * 0.453592 # [kg]

TAT_measured = np.array([13.8, 11.8, 9.2, 7.8, 6.8, 5.8])+ 273.15 # [Kelvin]

# Reference
# print("Reference")
# m = 9165 * 0.453592 + 95 + 92 + 74 + 66 + 61 + 75 + 78 + 86 + 68 + 4050 * 0.453592 # [kg]
# time_config1 = np.array([19*60+17, 21*60+37, 23*60+46, 26*60+4, 29*60+47, 32*60])
# IAS_config1 = (np.array([249, 221, 192, 163, 130, 118]) - 2) * 0.514444 # [m/s]
# pressure_alt = np.array([5010, 5020, 5020, 5030, 5020, 5110]) * 0.3048 # [m]
# AOA = np.array([1.7, 2.4, 3.6, 5.4, 8.7, 10.6])
# fuelflow_left = np.array([798, 673, 561, 463, 443, 474]) * 0.453592/60/60 # [kg/s]
# fuelflow_right = np.array([813, 682, 579, 484, 467, 499]) * 0.453592/60/60 # [kg/s]
# fuelused = np.array([360, 412, 447, 478, 532, 570]) * 0.453592 # [kg]

# TAT_measured = np.array([12.5, 10.5, 8.8, 7.2, 6, 5.2])+ 273.15 # [Kelvin]


AOA_rad = np.radians(AOA)

_lambda = -0.0065
T_0 = 288.15
gamma = 1.4

def p0overp(h_p):
    g0 = 9.80665
    R = 287.05
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
speedsound = np.sqrt(gamma * 287.05 * TAT_corrected)

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
V_TAS = Mach * speedsound
rho = 101325 / p0overp(pressure_alt) / 287.05 / TAT_corrected
Weight = (m - fuelused) * 9.80665

def Reynolds():
    # bulk_viscosity
    mu_0 = 1.7894e-5
    T_0 = 273.11
    _S = 110.56
    mu = mu_0 *( TAT_corrected / T_0 ) ** (1.5) * (T_0 + _S) / (TAT_corrected + _S)

    return rho * V_TAS * c / mu

Mach_min, Mach_max = np.min(Mach), np.max(Mach)
Reynolds_min, Reynolds_max = np.min(Reynolds()), np.max(Reynolds())

def CLCD_plot_stationary(plot = 'True'):

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

    cl_least = np.linspace(np.min(CL)-0.23, np.max(CL)+0.1,100)
    popt, pcov = spopt.curve_fit(func, CL, CD)
    CD0 = popt[0]
    oswald = 1/ popt[1] / np.pi / AR


    print('CD-CL ##################')
    print('a = ', popt[1])
    print('b or CD0 = ', CD0)
    print('oswald = ', oswald)


    ## Error derivation
    CD_leastsq = CD0 + CL**2 * popt[1]
    max_error = max(np.abs(CD-CD_leastsq))
    L2_error = np.sqrt( np.sum((CD-CD_leastsq)**2))
    error = max_error / L2_error * 100
    print('max error = ', max_error)
    print('L2 error =', L2_error)
    print('error =', error, '%')

    ## R2 analysis
    CD_avg = 1/np.shape(CD)[0] * np.sum(CD)
    SS_tot = np.sum((CD-CD_avg)**2)
    SS_res = np.sum((CD-CD_leastsq)**2)
    R2_CD = 1-SS_res/SS_tot
    print('R2 of CD = ', R2_CD)

    if plot == 'True':
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
        ax[0,0].plot(func(cl_least, *popt), cl_least , label = r'Least Square regression $C_D = aC_L^2 + b$')
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

def CLalpha_plot_stationary(plot='True'):

    CL = ( Weight - Thrust*np.sin(AOA_rad)) * 2 / (rho * V_TAS**2 * S)
    # CL = ( Weight ) * 2 / (rho * V_TAS**2 * S)

    ## Define text sizes for **SAVED** pictures (texpsize -- text export size)
    texpsize= [26,28,30]

    ## Input Arrays
    x = AOA
    y = CL

    ## best fit
    def func(AOA, a, b):
        return a*AOA + b

    popt, pcov = spopt.curve_fit(func, AOA, CL)
    CLalpha = popt[0]
    y_intercept = popt[1]
    aoa_least= np.linspace(np.min(AOA)-2, np.max(AOA)+0.1, 100)
    
    print('CD-alpha ################')
    print('a / CLalpha = ', CLalpha * 180 / np.pi)
    print('b = ', y_intercept)

    ## Error derivation
    CL_leastsq = y_intercept + AOA * CLalpha
    max_error = max(np.abs(CL-CL_leastsq))
    L2_error = np.sqrt( np.sum((CL-CL_leastsq)**2))
    error = max_error / L2_error * 100
    print('max error = ', max_error)
    print('L2 error =', L2_error)
    print('error =', error, '%')

    ## R2 analysis
    CL_avg = 1/np.shape(CL)[0] * np.sum(CL)
    SS_tot = np.sum((CL-CL_avg)**2)
    SS_res = np.sum((CL-CL_leastsq)**2)
    R2_CL = 1-SS_res/SS_tot
    print('R2 of CL = ', R2_CL)


    if plot == 'True':
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
        ax[0,0].plot(aoa_least, func(aoa_least, *popt), label = r'Least Square regression $C_L = a \alpha + b$')
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

    CLalpha_rad = CLalpha * 180 / np.pi
    return CLalpha_rad, y_intercept


# Theoretical plot
# def CL2CD_plot():
#     CL  = np.linspace(0,1.0,100)
#     CD = CLCD_plot_stationary()[0] + CL**2 / (np.pi * AR * CLCD_plot_stationary()[1])

#     ## Define text sizes for **SAVED** pictures (texpsize -- text export size)
#     texpsize= [26,28,30]

#     ## Input Arrays
#     x = np.linspace(1,10,10)
#     y = np.ones((10))*2

#     ## Graphing Parameters
#     SMALL_SIZE  = texpsize[0]
#     MEDIUM_SIZE = texpsize[1]
#     BIGGER_SIZE = texpsize[2]

#     plt.style.use('grayscale')
#     plt.rc('font', size=MEDIUM_SIZE, family='serif')    ## controls default text sizes
#     plt.rc('axes', titlesize=SMALL_SIZE)                ## fontsize of the axes title
#     plt.rc('axes', labelsize=SMALL_SIZE)                ## fontsize of the x and y labels
#     plt.rc('xtick', labelsize=SMALL_SIZE)               ## fontsize of the tick labels
#     plt.rc('ytick', labelsize=SMALL_SIZE)               ## fontsize of the tick labels
#     plt.rc('legend', fontsize=SMALL_SIZE)               ## legend fontsize
#     plt.rc('figure', titlesize=BIGGER_SIZE)             ## fontsize of the figure title
#     plt.rc('text', usetex=False)
#     matplotlib.rcParams['lines.linewidth']  = 1.5
#     matplotlib.rcParams['figure.facecolor'] = 'white'
#     matplotlib.rcParams['axes.facecolor']   = 'white'
#     matplotlib.rcParams["legend.fancybox"]  = False

#     ## Graph
#     fig, ax = plt.subplots(1,1,squeeze=False,figsize=(16,9))
#     ax[0,0].plot(CD, CL)
#     # ax[0,0].plot(x+x, y+y, label="test2", linestyle="dashed")
#     #ax[0,0].loglog(x, y, marker = "s", color='black', markerfacecolor='none', markeredgewidth=2, markersize=6, label="test")
#     ax[0,0].set_ylabel(r"$C_D$")          ## String is treatable as latex code
#     ax[0,0].set_xlabel(r"$C_L$")
#     #ax[0,0].set_xlim(0,x[-1])
#     ax[0,0].grid(True,which="major",color="#999999")
#     ax[0,0].grid(True,which="minor",color="#DDDDDD",ls="--")
#     ax[0,0].minorticks_on()
#     ax[0,0].tick_params(which='major', length=10, width=2, direction='inout')
#     ax[0,0].tick_params(which='minor', length=5, width=2, direction='in')
#     # ax[0,0].legend(loc=0, framealpha=1.0).get_frame().set_edgecolor('k')
#     fig.savefig("CLCD-theoretical.png", bbox_inches='tight')                                    ## Insert save destination

#     ## If you want to see the figure, else disable last two lines.
#     fig.tight_layout()
#     plt.show()

#     return 0

# CLCD_plot_stationary()
# CLalpha_plot_stationary()
# CL2CD_plot()

# variables to be exported
CD0, oswald = CLCD_plot_stationary('False')
CLalpha, y_intercept = CLalpha_plot_stationary('False')
CL = ( Weight - Thrust*np.sin(AOA_rad)) * 2 / (rho * V_TAS**2 * S)





