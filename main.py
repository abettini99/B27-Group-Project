"""
Institution:    TU Delft
Authors:        B27
Date:           06-03-2020
Dynamic Analysis of Cessna Citation II
"""

# ==============================================================================================
# Import Libraries
# ==============================================================================================
from math import pi, pow, sin, cos, radians, degrees    # provides access to the mathematical functions defined by the C standard
import pandas as pd                                     # package for improved data analysis through DataFrames, etc...
import numpy as np                                      # fundamental package for scientific computing
import control.matlab as ml                             # import module to emulate functionality of MATLAB
import control as ctl                                   # import package for analysis and design of feedback control systems
import matplotlib.pyplot as plt                         # package to create visualisations
from scipy.io import loadmat                            # loadmat imports a .mat file
from numpy.linalg import inv, eig                       # inv computes the inverse of a matrix; eig computes eigenvalues of matrix
import matplotlib                                       # comprehensive library for creating static, animated, and interactive visualizations in Python.

# ==============================================================================================
# Function Definitions
# ==============================================================================================
def importdata(filename):
    """
        This function imports the data from the flight tests
        :filename: relative path of the .mat-file from the flight test
        :return: dataframe with each variable in one column
    """
    mat     = loadmat(filename)          # load data from .mat-file - returns dictionary with variable names as keys and loaded matrices as values
    mdata   = mat['flightdata']          # access first level variable in .mat-file
    mdtype  = mdata.dtype                # dtypes of structures are unsized objects
    ndata   = {n: mdata[n][0,0]['data'][0,0] for n in mdtype.names}                       # redefine data as dictionary with names from dtypes
    data    = pd.DataFrame(data=None, index=None, columns=None, dtype=None, copy=False)   # initialise empty dataframe

    for key, values in ndata.items():    # iterate over keys and
        data[key] = values.flatten()     # add new column with variable as key and values; input to dataframe must be 1D such that 2D arrays must be flattened
    return data

def manouvre(data, flightmanouvre):
    """
        This function slices the dataframe into a smaller dataframe for each flight manouvre with the corresponding start and stop time
        :flightmanouvre: name of flightmanouvre (phugoid, shortperiodoscillation, heavilydampedmotion, spiral or dutchroll)
        :return: sliced dataframe with each variable in one column
    """
    if flightmanouvre == "clcd":
        time_start  = 992
        time_stop   = 1740
        data        = data[(data['time'] >= time_start) & (data['time'] <= time_stop)]
        data        = data[data.measurement_running != 0]
        return data

    if flightmanouvre == "elevatortrim":
        time_start  = 1800
        time_stop   = 2340
        data        = data[(data['time'] >= time_start) & (data['time'] <= time_stop)]
        data        = data[data.measurement_running != 0]
        return data

    if flightmanouvre == "cgshift":
        time_start  = 2340
        time_stop   = 2640
        data        = data[(data['time'] >= time_start) & (data['time'] <= time_stop)]
        data        = data[data.measurement_running != 0]
        return data

    if flightmanouvre == "phugoid":
        time_start  = 2675
        time_stop   = 2810
        data        = data[(data['time'] >= time_start) & (data['time'] <= time_stop)]
        return data

    if flightmanouvre == "shortperiod":
        time_start  = 2636
        time_stop   = 2670
        data        = data[(data['time'] >= time_start) & (data['time'] <= time_stop)]
        return data

    if flightmanouvre == "dutchroll":
        time_start  = 3020
        time_stop   = 3041.5
        data        = data[(data['time'] >= time_start) & (data['time'] <= time_stop)]
        return data
    if flightmanouvre == "dutchrollYD":
        time_start  = 3091
        time_stop   = 3103
        data        = data[(data['time'] >= time_start) & (data['time'] <= time_stop)]
        return data

    if flightmanouvre == "aperroll":
        time_start  = 2880
        time_stop   = 2880 + 150
        data        = data[(data['time'] >= time_start) & (data['time'] <= time_stop)]
        return data

    if flightmanouvre == "spiral":
        time_start  = 3300-10
        time_stop   = 3507
        data        = data[(data['time'] >= time_start) & (data['time'] <= time_stop)]
        return data

def cgshift1(mi, cgi, mf, momentchange):
    """
        Function description
        :mi:
        :cgi:
        :mf:
        :momentchange:
        :return:
    """
    return (mi * cgi+ momentchange) / mf

def cgshift2(m, cgi, momentchange):
    """
        Function description
        :m:
        :cgi:
        :momentchange:
        :return:
    """
    return cgi + momentchange / m

def cgmac(cgdatum, lemac, mac):
    """
        Function description
        :cgdatum:
        :lemac: leading edge position of mean aerodynamic chord
        :mac: mean aerodynamic chord
        :return:
    """
    return (cgdatum - lemac) / mac * 100

def fuelinterpolation(fuel, fuelmoment):
    """
        Function description
        :fuel:
        :fuelmoment:
        :return:
    """
    for i in range(len(fuelmoment)):
        if fuel<= fuelmoment.iat[i,0]:
            momentcg = fuelmoment.iat[i,1] - ((fuelmoment.iat[i,1] - fuelmoment.iat[i-1,1])/(fuelmoment.iat[i,0]-fuelmoment.iat[i-1,0]) )* (fuelmoment.iat[i,0] - fuel)
            break
    return momentcg

def trapezoidal(j, time, ff_le, ff_re):
    """
        Function description
        :j:
        :time:
        :ff_le:
        :ff_re:
        :return:
    """
    dt = time[j+1] -time[j]
    f1 = ff_le[j] / 3600
    f2 = ff_le[j+1] / 3600
    f3 = ff_re[j] / 3600
    f4 = ff_re[j+1] / 3600
    I1 = (f1+f2) * dt / 2
    I2 = (f3+f4) * dt / 2
    return I1+I2

def fuelused(time, ff_le, ff_re):
    """
        Function description
        :time:
        :ff_le:
        :ff_re:
        :return:
    """
    fuelUsed = {}
    integral = 0
    for t in time[:len(time)-1]:
        integral =integral + trapezoidal(time[time == t].index[0], time, ff_le, ff_re)
        fuelUsed[t] = integral
    return fuelUsed

def cgtime(t, fuel_i, ZFM, CG_ZFM, fuelused, fuel_moment):
    """
        Function description
        :t:
        :fuel_i:
        :ZFM: zero fuel mass [UNIT]
        :CG_ZFM:
        :fuelused:
        :fuel_moment:
        :return:
    """
    fuel_mass = fuel_i - fuelused[t]
    fuel_mom = fuelinterpolation(fuel_mass, fuel_moment)*100
    tot_mass = ZFM+ fuel_mass
    cg = cgshift1(ZFM, CG_ZFM, tot_mass, fuel_mom)
    return cg, tot_mass

def ktstoms(velocity):
    """
        Function converts velocity given in knots to metre per second
        :velocity: input velocity [kts]
        :return: output velocity [m/s]
    """
    return velocity * 0.5144444444

def fttom(altitude):
    """
        Function converts altitude given in feet [ft] to metre [m]
        :velocity: input altitude in [ft]
        :return: outputs altitude in [m]
    """
    return altitude * 0.3048

def lbstokg(mass):
    """
        Function converts mass given in pounds [lbs] to kilogram [kg]
        :velocity: input mass in [lbs]
        :return: outputs mass in [kg]
    """
    return mass * 0.45359237

# ==============================================================================================
# Import data from Matlab files
# ==============================================================================================
# data = importdata('referencedata.mat')                                                     # initialise reference data from matlab file
data = importdata('flightdata.mat')                                                          # initialise flight data from matlab file

# ==============================================================================================
# Stationary measurements
# ==============================================================================================
clcd    = manouvre(data, 'clcd')                                                             # sliced data for the 6 CL-CD measurement series
etrim   = manouvre(data, 'elevatortrim')                                                     # sliced data for the 7 e-trim measurement series
cgshift = manouvre(data, 'cgshift')                                                          # sliced data for the 2 cg-shift measurement series

# ==============================================================================================
# centre of gravity calculations
# ==============================================================================================
momentfuel = pd.read_excel('FuelCG.xlsx', header=None, sheet_name='Sheet1')

LEMAC   = 261.56*0.0254                                                                  # [m] leading edge position of mean aerodynamic chord
BEM     = 9165                                                                           # [lbs] basic empty mass, taken from weight measurements
cgBEM   = 291.647954                                                                     # [in] centre of gravity position of BEM

xseats  = [131,131,214,214,251,251,288,288,170]                                          # [in] x-position of each seat for passenger / pilot
Mseats  = [90*2.20462,102*2.20462,83*2.20462,94*2.20462,84*2.20462, \
           74*2.20462,79*2.20462,103*2.20462,80*2.20462]                                 # [lbs] weight of each passenger on each seat

xbag    = [74, 321, 338]                                                                 # [in]
Mbag    = [0, 0, 0]                                                                      # [lbs]

mPL     = sum(Mseats) + sum(Mbag)                                                        # [lbs] payload mass

momentPL = sum([xseats[i]*Mseats[i] for i in range(len(xseats))]) \
           + sum([xbag[j]*Mbag[j] for j in range(len(xbag))])

ZFM = BEM + mPL                                                                          # [lbs] zero fuel mass
cgZFM  = cgshift1(BEM, cgBEM, ZFM, momentPL)                                             # [] centre of gravity position of zero fuel mass

mf_init = 4100                                                                           # [lbs] initial fuel mass
momentmf_init = fuelinterpolation(mf_init, momentfuel)*100

mRW    = ZFM + mf_init                                                                   # [lbs] ramp weight
cgmRW  = cgshift1(ZFM, cgZFM, mRW, momentfuel)                                           # [] centre of gravity of ramp weight

cg     = []                                                                              # initialise empty list for centre of gravity
temp   = np.zeros((len(data.time)-1, 2))                                                 # initialise empty numpy array for time and mass
fused  = fuelused(data.time, data.rh_engine_FMF, data.lh_engine_FMF)
i = 0
for t in data.time[:len(data.time) - 1]:
    tempcg, tempm = cgtime(t, mf_init, ZFM, cgZFM, fused, momentfuel)
    cg.append(tempcg)                                                                    # append invidual cg position to cg list
    temp[i][0], temp[i][1] = t, lbstokg(tempm)                                           # population of numpy array with time and mass
    i += 1

cg     = np.array(cg)                                                                    # convert list to numpy array
df     = pd.DataFrame(temp, columns=['time', 'mass'])                                    # [s, kg] dataframe with specific mass at time t
# ==============================================================================================
# Eigenmotion analysis
# ==============================================================================================
f = open('eigenvalues.txt', 'w+')                                                            # create .txt-file where EV's are written

for motion in ['phugoid', 'shortperiod', 'aperroll', 'dutchroll', 'dutchrollYD', 'spiral']:
    # ==============================================================================================
    # Import data from Matlab files
    # ==============================================================================================
    # data = importdata('referencedata.mat')                                                   # initialise reference data from matlab file
    data   = importdata('flightdata.mat')                                                    # initialise flight data from matlab file
    data   = manouvre(data, motion)                                                          # sliced data array for phugoid motion
    m      = df
    m      = manouvre(m, motion)

    # ==============================================================================================
    # Transform coordinate system from body to stability axis frame
    # ==============================================================================================
    alpha0 = radians(data.vane_AOA.iloc[0])                                                  # [rad] angle of attack in the stationary flight condition
    theta0 = radians(data.Ahrs1_Pitch.iloc[0])                                               # [rad] pitch angle in the stationary flight condition

    data['vane_AOA']    = data['vane_AOA'] - degrees(alpha0)                                 # Transform angle of attack from body to stability axis frame
    data['Ahrs1_Pitch'] = data['Ahrs1_Pitch'] - degrees(theta0)                              # Transform angle of pitch from body to stability axis frame

    # ==============================================================================================
    # Parameter definition; copied from Cit_par.py
    # ==============================================================================================
    hp0    = fttom(data.Dadc1_alt.iloc[0])                                                   # [m] pressure altitude in the stationary flight condition
    V0     = ktstoms(data.Dadc1_tas.iloc[0])                                                 # [m/s] true airspeed in the stationary flight condition

    m      = m.mass.iloc[0]                                                                  # [kg] takeoff weight of Cessna Citation II

    e      = 0.5662173213929709                                                              # [-] Oswald factor
    # e      = 0.8                                                                             # [-] Oswald factor reference data
    CD0    = 0.0255504796766271                                                              # [-] Zero lift drag coefficient
    # CD0    = 0.04                                                                            # [-] Zero lift drag coefficient reference data
    CLa    = 4.557979876353312                                                               # [1/rad] Slope of CL-alpha curve
    # CLa    = 5.084                                                                           # [1/rad] Slope of CL-alpha curve reference data

    Cma    = -0.582128                                                                       # [-] longitudinal stabilty
    Cmde   = -1.21076                                                                        # [-] elevator effectiveness

    S      = 30.00                                                                           # [m^2] wing area
    Sh     = 0.2 * S                                                                         # [m^2] stabiliser area
    Sh_S   = Sh / S                                                                          # [-]
    lh     = 0.71 * 5.968                                                                    # [m] tail length
    c      = 2.0569                                                                          # [m] mean aerodynamic cord
    lh_c   = lh / c                                                                          # [-]
    b      = 15.911                                                                          # [m] wing span
    bh     = 5.791                                                                           # [m] stabilser span
    AR     = b ** 2 / S                                                                      # [-] wing aspect ratio
    Ah     = bh ** 2 / Sh                                                                    # [-] stabilser aspect ratio
    Vh_V   = 1                                                                               # [-]
    ih     = -2 * pi / 180                                                                   # [rad] stabiliser angle of incidence

    rho0   = 1.2250                                                                          # [kg/m^3] air density at sea level
    LAMBDA = -0.0065                                                                         # [K/m] temperature gradient in ISA
    Temp0  = 288.15                                                                          # [K] temperature at sea level in ISA
    R      = 287.05                                                                          # [m^2/s^2 K] specific gas constant
    g      = 9.81                                                                            # [m/s^2] gravitational acceleration
    rho    = rho0 * pow( ((1+(LAMBDA * hp0 / Temp0))), (-((g / (LAMBDA*R)) + 1)))            # [kg/m^3] density at altitude h
    W      = m * g                                                                           # [N] aircraft weight

    muc    = m / (rho * S * c)
    mub    = m / (rho * S * b)
    KX2    = 0.019
    KY2    = 1.3925
    KZ2    = 0.042
    KXZ    = 0.002

    Cmac   = 0                                                                               # [-] Moment coefficient about the aerodynamic centre
    CNwa   = CLa                                                                             # [-] Wing normal force slope
    CNha   = 2 * pi * Ah / (Ah + 2)                                                          # [-] Stabiliser normal force slope
    depsda = 4 / (AR + 2)                                                                    # [-] Downw sh gradient

    CL     = 2 * W / (rho * V0 ** 2 * S)                                                     # [-] Lift coefficient
    CD     = CD0 + (CLa * alpha0) ** 2 / (pi * AR * e)                                       # [-] Drag coefficient

    CX0    = W * sin(theta0) / (0.5 * rho * V0 ** 2 * S)
    CXu    = -0.095
    CXa    = +0.47966
    CXadot = +0.08330
    CXq    = -0.28170
    CXde   = -0.03728

    CZ0    = -W * cos(theta0) / (0.5 * rho * V0 ** 2 * S)
    CZu    = -0.37616
    CZa    = -5.74340
    CZadot = -0.00350
    CZq    = -5.66290
    CZde   = -0.69612

    Cmu    = +0.06990
    Cmadot = +0.17800
    Cmq    = -8.79415

    CYb    = -0.7500
    CYbdot = +0
    CYp    = -0.0304
    CYr    = +0.8495
    CYda   = -0.0400
    CYdr   = +0.2300

    Clb    = -0.10260
    Clp    = -0.71085
    Clr    = +0.23760
    Clda   = -0.23088
    Cldr   = +0.03440

    Cnb    =  +0.1348
    Cnbdot =  +0
    Cnp    =  -0.0602
    Cnr    =  -0.2061
    Cnda   =  -0.0120
    Cndr   =  -0.0939

    # ==============================================================================================
    # Declaration of matrices and column vectors
    # ==============================================================================================
    A       = np.zeros((8,8))                                                                # Declaration of matrix A with dimensions [8 x 8] for system of equations
    B       = np.zeros((8,8))                                                                # Declaration of matrix B with dimensions [8 x 4] for system of equations
    C       = np.zeros((8,4))                                                                # Declaration of matrix C with dimensions [8 x 4] for system of equations
    A_temp  = np.zeros((8,8))                                                                # Declaration of temporary matrix A with dimensions [8 x 8] for system of equations
    B_temp  = np.zeros((8,8))                                                                # Declaration of temporary matrix B with dimensions [8 x 4] for system of equations
    C_temp  = np.zeros((8,4))                                                                # Declaration of temporary matrix C with dimensions [8 x 4] for system of equations
    As      = np.zeros((4,4))                                                                # Declaration of matrix As with dimensions [4 x 4] for symmetric EOM
    Aa      = np.zeros((4,4))                                                                # Declaration of matrix Aa with dimensions [4 x 4] for asymmetric EOM
    Bs      = np.zeros((4,4))                                                                # Declaration of matrix Bs with dimensions [4 x 2] for symmetric EOM
    Ba      = np.zeros((4,4))                                                                # Declaration of matrix Ba with dimensions [4 x 2] for asymmetric EOM
    Cs      = np.zeros((4,2))                                                                # Declaration of matrix Cs with dimensions [4 x 2] for symmetric EOM
    Ca      = np.zeros((4,2))                                                                # Declaration of matrix Ca with dimensions [4 x 2] for asymmetric EOM

    # ==============================================================================================
    # Population of symmetric EOM matrices with variables for state-space representation
    # ==============================================================================================
    As[0,0] = - 2 * muc * c / V0**2

    As[1,1] = (CZadot - 2 * muc) * c / V0

    As[2,2] = -c / V0

    As[3,1] = Cmadot * c / V0
    As[3,3] = -2 * muc * KY2 * (c / V0)**2

    Bs[0,0] = -CXu / V0
    Bs[0,1] = -CXa
    Bs[0,2] = -CZ0
    Bs[0,3] = -CXq * c / V0

    Bs[1,0] = -CZu / V0
    Bs[1,1] = -CZa
    Bs[1,2] = +CX0
    Bs[1,3] = -(CZq + 2 * muc) * c / V0

    Bs[2,3] = -c / V0

    Bs[3,0] = -Cmu / V0
    Bs[3,1] = -Cma
    Bs[3,3] = -Cmq * c / V0

    Cs[0,0] = -CXde

    Cs[1,0] = -CZde

    Cs[3,0] = -Cmde

    # ==============================================================================================
    # Population of asymmetric EOM matrices with variables for state-space representation
    # ==============================================================================================
    Aa[0,0] = (CYbdot - 2 * mub) * b / V0

    Aa[1,1] = -0.5 * b / V0

    Aa[2,2] = -4 * mub * KX2 * b**2 / (2 * V0**2)
    Aa[2,3] = 4 * mub * KXZ * b**2 / (2 * V0**2)

    Aa[3,0] = Cnbdot * b / V0
    Aa[3,2] = 4 * mub * KXZ * b**2 / (2 * V0**2)
    Aa[3,3] = -4 * mub * KZ2 * b**2 / (2 * V0**2)

    Ba[0,0] = -CYb
    Ba[0,1] = -CL
    Ba[0,2] = -CYp * b / (2 * V0)
    Ba[0,3] = -(CYr - 4 * mub) * b / (2 * V0)

    Ba[1,2] = -b / (2 * V0)

    Ba[2,0] = -Clb
    Ba[2,2] = -Clp * b / (2 * V0)
    Ba[2,3] = -Clr * b / (2 * V0)

    Ba[3,0] = -Cnb
    Ba[3,2] = -Cnp * b / (2 * V0)
    Ba[3,3] = -Cnr * b / (2 * V0)

    Ca[0,0] = -CYda
    Ca[0,1] = -CYdr

    Ca[2,0] = -Clda
    Ca[2,1] = -Cldr

    Ca[3,0] = -Cnda
    Ca[3,1] = -Cndr

    # ==============================================================================================
    # Population of matrices As and Aa
    # ==============================================================================================
    A_temp[0:4, 0:4] = As
    A_temp[4:8, 4:8] = Aa
    B_temp[0:4, 0:4] = Bs
    B_temp[4:8, 4:8] = Ba
    C_temp[0:4, 0:2] = Cs
    C_temp[4:8, 2:4] = Ca

    As = np.dot(inv(As), Bs)
    Bs = np.dot(inv(As), Cs)

    Aa = np.dot(inv(Aa), Ba)
    Ba = np.dot(inv(Aa), Ca)

    # Output of state-space representation should be equal to the relevant aircraft states
    # --> matrix C is the identity matrix and D is a zero array
    C = np.identity(4)
    D = np.zeros((4,2))

    # ==============================================================================================
    # Calculates responses to symmetric eigenmotions from state-space system
    # ==============================================================================================
    if motion in ['phugoid', 'shortperiod']:
        syss = ctl.StateSpace(As, Bs, C, D)                                                  # create state-space system for symmetric eigenmotions
        evals, evacs = eig(As)                                                               # compute eigenvalues and eigenvectors

        for i in range(0, len(evals)):                                                       # write eigenvalues to textfile
            f = open('eigenvalues.txt', 'a+')                                                # append lines to existing .txt-file
            f.write("{}, lambda{}: {} \n".format(motion, (i+1), evals[i]))                   # write eigenvalues

        tstop = data.time.iloc[-1] - data.time.iloc[0]                                       # normalise final time value for manouvre
        dt  = np.arange(0, tstop + 0.1, 0.1)                                                 # create time vector with 0.1s step size

        units = ['[m/s]', '[rad]', '[rad]', '[rad/s]']                                       # list with units of columns for plotting
        u = [np.radians(data.delta_e), np.zeros(len(data.index))]                            # [rad] input array given input at each time for [de, dt]
        columns = [r'V_{TAS}', r'\alpha', r'\theta', r'q']                                   # names of invidiual columns for DataFrame
        eigenmotion = []                                                                     # initialise empty list

        flightdata = [np.radians(data.vane_AOA), np.radians(data.Ahrs1_Pitch), np.radians(data.Ahrs1_bPitchRate)]

        t, y, x = ctl.forced_response(syss, dt, U=u)                                         # calculate forced response
        df2 = pd.DataFrame(np.transpose(y), columns=columns)                                 # convert forced response to DataFrame
        eigenmotion.append(df2)                                                              # append DataFrame to individual list
        eigenmotion = pd.concat(eigenmotion, axis=1)                                         # concatenate list into panda dataframe along axis 1

        if motion == 'phugoid':
            fig1, ax1 = plt.subplots(4,1, squeeze=False, figsize=(16,9))                     # initialise figure with 4 rows and 1 column
            for i in range(0,3):
                ax1[i,0].plot(t, eigenmotion.iloc[:,i+1], 'C1', label='Numerical Model')     # plot each variable from output vector
                ax1[i,0].plot(t, flightdata[i], c='k', label='Experimental Data')            # plot each variable from test flight data
                ax1[i,0].set_xlabel('$t$ [s]')                                               # set label of x-axis
                ax1[i,0].set_xlim(xmin=0)                                                    # set xmin at 0
                ax1[i,0].set_ylabel('${}$ {}'.format(eigenmotion.columns[i+1], units[i+1]))  # set label of y-axis
                ax1[i,0].minorticks_on()                                                     # set minor ticks
                ax1[i,0].grid(which='major', linestyle='-', linewidth='0.5', color='black')  # customise major grid
                ax1[i,0].grid(which='minor', linestyle=':', linewidth='0.5', color='grey')   # customise minor grid
                ax1[i,0].legend(loc=0, framealpha=1.0).get_frame().set_edgecolor('k')        # set legend for subplot

            ax1[3,0].plot(t, u[0], c='k', label='Elevator Deflection')                       # plot input variable
            ax1[3,0].set_xlabel('$t$ [s]')                                                   # set label of x-axis
            ax1[3,0].set_xlim(xmin=0)                                                        # set xmin at 0
            ax1[3,0].set_ylabel('$\delta_e$ [rad]')                                          # set label of y-axis
            ax1[3,0].minorticks_on()                                                         # set minor ticks
            ax1[3,0].grid(which='major', linestyle='-', linewidth='0.5', color='black')      # customise major grid
            ax1[3,0].grid(which='minor', linestyle=':', linewidth='0.5', color='grey')       # customise minor grid
            ax1[3,0].legend(loc=0, framealpha=1.0).get_frame().set_edgecolor('k')            # set legend

            fig1.tight_layout(pad=1.0)                                                       # increase spacing between subplots
            # fig1.suptitle('Phugoid')                                                         # set title of figure
            fig1.savefig('images/phugoid.png', dpi=300, bbox_inches='tight')                 # save figure

        elif motion == 'shortperiod':
            fig1, ax1 = plt.subplots(4,1, squeeze=False, figsize=(16,9))                     # initialise figure with 4 rows and 1 column
            for i in range(0,3):
                ax1[i,0].plot(t, eigenmotion.iloc[:,i+1], 'C21', label='Numerical Model')    # plot each variable from output vector
                ax1[i,0].plot(t, flightdata[i], c='k', label='Experimental Data')            # plot each variable from test flight data
                ax1[i,0].set_xlabel('$t$ [s]')                                               # set label of x-axis
                ax1[i,0].set_xlim(xmin=0)                                                    # set xmin at 0
                ax1[i,0].set_ylabel('${}$ {}'.format(eigenmotion.columns[i+1], units[i+1]))  # set label of y-axis
                ax1[i,0].minorticks_on()                                                     # set minor ticks
                ax1[i,0].grid(which='major', linestyle='-', linewidth='0.5', color='black')  # customise major grid
                ax1[i,0].grid(which='minor', linestyle=':', linewidth='0.5', color='grey')   # customise minor grid
                ax1[i,0].legend(loc=0, framealpha=1.0).get_frame().set_edgecolor('k')        # set legend for subplot

            ax1[3,0].plot(t, u[0], c='k', label='Elevator Deflection')                       # plot input variable
            ax1[3,0].set_xlabel('$t$ [s]')                                                   # set label of x-axis
            ax1[3,0].set_xlim(xmin=0)                                                        # set xmin at 0
            ax1[3,0].set_ylabel('$\delta_e$ [rad]')                                          # set label of y-axis
            ax1[3,0].minorticks_on()                                                         # set minor ticks
            ax1[3,0].grid(which='major', linestyle='-', linewidth='0.5', color='black')      # customise major grid
            ax1[3,0].grid(which='minor', linestyle=':', linewidth='0.5', color='grey')       # customise minor grid
            ax1[3,0].legend(loc=0, framealpha=1.0).get_frame().set_edgecolor('k')            # set legend

            fig1.tight_layout(pad=1.0)                                                       # increase spacing between subplots
            # fig1.suptitle('Short Period Oscillation')                                        # set title of figure
            fig1.savefig('images/shortperiod.png', dpi=300, bbox_inches='tight')             # save figure

    # ==============================================================================================
    # Calculates responses to asymmetric eigenmotions from state-space system
    # ==============================================================================================
    if motion in ['aperroll', 'dutchroll', 'dutchrollYD', 'spiral']:
        sysa = ctl.StateSpace(Aa, Ba, C, D)                                                  # create state-space system for symmetric eigenmotions
        evals, evacs = eig(Aa)                                                               # compute eigenvalues and eigenvectors

        for i in range(0, len(evals)):                                                       # write eigenvalues to textfile
            f = open('eigenvalues.txt', 'a+')                                                # append lines to existing .txt-file
            f.write("{}, lambda{}: {} \n".format(motion, (i+1), evals[i]))                   # write eigenvalues

        tstop = data.time.iloc[-1] - data.time.iloc[0]                                       # normalise final time value for manouvre
        dt  = np.arange(0, tstop + 0.1, 0.1)                                                 # create time vector with 0.1s step size

        units = ['[rad]', '[rad]', '[rad/s]', '[rad/s]']                                     # list with units of columns for plotting
        u = [np.radians(data.delta_a), np.radians(data.delta_r)]                             # [rad] input array given input at each time for [da, dr]
        columns = [r'\beta', r'\phi', r'p', r'r']                                            # names of invidiual columns for DataFrame
        eigenmotion = []                                                                     # initialise empty list

        flightdata = [np.radians(data.Ahrs1_Roll), np.radians(data.Ahrs1_bRollRate), np.radians(data.Ahrs1_bYawRate)]

        t, y, x = ctl.forced_response(sysa, dt, U=u)                                         # calculate forced response
        df2 = pd.DataFrame(np.transpose(y), columns=columns)                                 # convert forced response to DataFrame
        eigenmotion.append(df2)                                                              # append DataFrame to individual list
        eigenmotion = pd.concat(eigenmotion, axis=1)                                         # concatenate list into panda dataframe along axis 1

        if motion == 'aperroll':
            fig1, ax1 = plt.subplots(4,1, squeeze=False, figsize=(16,9))                     # initialise figure with 4 rows and 1 column
            for i in range(0,3):
                ax1[i,0].plot(t, eigenmotion.iloc[:,i+1], 'C1', label='Numerical Model')     # plot each variable from output vector
                ax1[i,0].plot(t, flightdata[i], c='k', label='Experimental Data')            # plot each variable from test flight data
                ax1[i,0].set_xlabel('$t$ [s]')                                               # set label of x-axis
                ax1[i,0].set_xlim(xmin=0)                                                    # set xmin at 0
                ax1[i,0].set_ylabel('${}$ {}'.format(eigenmotion.columns[i+1], units[i+1]))  # set label of y-axis
                ax1[i,0].minorticks_on()                                                     # set minor ticks
                ax1[i,0].grid(which='major', linestyle='-', linewidth='0.5', color='black')  # customise major grid
                ax1[i,0].grid(which='minor', linestyle=':', linewidth='0.5', color='grey')   # customise minor grid
                ax1[i,0].legend(loc=0, framealpha=1.0).get_frame().set_edgecolor('k')        # set legend for subplot

            ax1[3,0].plot(t, u[0], c='k', linestyle='--', label='Aileron Deflection')        # plot input variable delta_a
            ax1[3,0].plot(t, u[1], c='k', linestyle='-',label='Rudder Deflection')           # plot input variable delta_r
            ax1[3,0].set_xlabel('$t$ [s]')                                                   # set label of x-axis
            ax1[3,0].set_xlim(xmin=0)                                                        # set xmin at 0
            ax1[3,0].set_ylabel('$\delta_a, \delta_r$ [rad]')                                # set label of y-axis
            ax1[3,0].minorticks_on()                                                         # set minor ticks
            ax1[3,0].grid(which='major', linestyle='-', linewidth='0.5', color='black')      # customise major grid
            ax1[3,0].grid(which='minor', linestyle=':', linewidth='0.5', color='grey')       # customise minor grid
            ax1[3,0].legend(loc=0, framealpha=1.0).get_frame().set_edgecolor('k')            # set legend

            fig1.tight_layout(pad=1.0)                                                       # increase spacing between subplots
            # fig1.suptitle('Aperiodic Roll')                                                  # set title of figure
            fig1.savefig('images/aperiodicroll.png', dpi=300, bbox_inches='tight')           # save figure

        if motion == 'dutchroll':
            fig1, ax1 = plt.subplots(4,1, squeeze=False, figsize=(16,9))                     # initialise figure with 4 rows and 1 column
            for i in range(0,3):
                ax1[i,0].plot(t, eigenmotion.iloc[:,i+1], 'C1', label='Numerical Model')     # plot each variable from output vector
                ax1[i,0].plot(t, flightdata[i], c='k', label='Experimental Data')            # plot each variable from test flight data
                ax1[i,0].set_xlabel('$t$ [s]')                                               # set label of x-axis
                ax1[i,0].set_xlim(xmin=0)                                                    # set xmin at 0
                ax1[i,0].set_ylabel('${}$ {}'.format(eigenmotion.columns[i+1], units[i+1]))  # set label of y-axis
                ax1[i,0].minorticks_on()                                                     # set minor ticks
                ax1[i,0].grid(which='major', linestyle='-', linewidth='0.5', color='black')  # customise major grid
                ax1[i,0].grid(which='minor', linestyle=':', linewidth='0.5', color='grey')   # customise minor grid
                ax1[i,0].legend(loc=0, framealpha=1.0).get_frame().set_edgecolor('k')        # set legend for subplot

            ax1[3,0].plot(t, u[0], c='k', linestyle='--', label='Aileron Deflection')        # plot input variable delta_a
            ax1[3,0].plot(t, u[1], c='k', linestyle='-',label='Rudder Deflection')           # plot input variable delta_r
            ax1[3,0].set_xlabel('$t$ [s]')                                                   # set label of x-axis
            ax1[3,0].set_xlim(xmin=0)                                                        # set xmin at 0
            ax1[3,0].set_ylabel('$\delta_a, \delta_r$ [rad]')                                # set label of y-axis
            ax1[3,0].minorticks_on()                                                         # set minor ticks
            ax1[3,0].grid(which='major', linestyle='-', linewidth='0.5', color='black')      # customise major grid
            ax1[3,0].grid(which='minor', linestyle=':', linewidth='0.5', color='grey')       # customise minor grid
            ax1[3,0].legend(loc=0, framealpha=1.0).get_frame().set_edgecolor('k')            # set legend

            fig1.tight_layout(pad=1.0)                                                       # increase spacing between subplots
            # fig1.suptitle('Dutch Roll without Yaw Damper')                                   # set title of figure
            fig1.savefig('images/dutchroll.png', dpi=300, bbox_inches='tight')               # save figure

        if motion == 'dutchrollYD':
            fig1, ax1 = plt.subplots(4,1, squeeze=False, figsize=(16,9))                     # initialise figure with 4 rows and 1 column
            for i in range(0,3):
                ax1[i,0].plot(t, eigenmotion.iloc[:,i+1], 'C1', label='Numerical Model')     # plot each variable from output vector
                ax1[i,0].plot(t, flightdata[i], c='k', label='Experimental Data')            # plot each variable from test flight data
                ax1[i,0].set_xlabel('$t$ [s]')                                               # set label of x-axis
                ax1[i,0].set_xlim(xmin=0)                                                    # set xmin at 0
                ax1[i,0].set_ylabel('${}$ {}'.format(eigenmotion.columns[i+1], units[i+1]))  # set label of y-axis
                ax1[i,0].minorticks_on()                                                     # set minor ticks
                ax1[i,0].grid(which='major', linestyle='-', linewidth='0.5', color='black')  # customise major grid
                ax1[i,0].grid(which='minor', linestyle=':', linewidth='0.5', color='grey')   # customise minor grid
                ax1[i,0].legend(loc=0, framealpha=1.0).get_frame().set_edgecolor('k')        # set legend for subplot

            ax1[3,0].plot(t, u[0], c='k', linestyle='--', label='Aileron Deflection')        # plot input variable delta_a
            ax1[3,0].plot(t, u[1], c='k', linestyle='-',label='Rudder Deflection')           # plot input variable delta_r
            ax1[3,0].set_xlabel('$t$ [s]')                                                   # set label of y-axis
            ax1[3,0].set_xlim(xmin=0)                                                        # set xmin at 0
            ax1[3,0].set_ylabel('$\delta_a, \delta_r$ [rad]')                                # set label of y-axis
            ax1[3,0].minorticks_on()                                                         # set minor ticks
            ax1[3,0].grid(which='major', linestyle='-', linewidth='0.5', color='black')      # customise major grid
            ax1[3,0].grid(which='minor', linestyle=':', linewidth='0.5', color='grey')       # customise minor grid
            ax1[3,0].legend(loc=0, framealpha=1.0).get_frame().set_edgecolor('k')            # set legend

            fig1.tight_layout(pad=1.0)                                                       # increase spacing between subplots
            # fig1.suptitle('Dutch Roll with Yaw Damper')                                      # set title of figure
            fig1.savefig('images/dutchrollYD.png', dpi=300, bbox_inches='tight')             # save figure

        if motion == 'spiral':
            fig1, ax1 = plt.subplots(4,1, squeeze=False, figsize=(16,9))                     # initialise figure with 4 rows and 1 column
            for i in range(0,3):
                ax1[i,0].plot(t, eigenmotion.iloc[:,i+1], 'C1', label='Numerical Model')     # plot each variable from output vector
                ax1[i,0].plot(t, flightdata[i], c='k', label='Experimental Data')            # plot each variable from test flight data
                ax1[i,0].set_xlabel('$t$ [s]')                                               # set label of x-axis
                ax1[i,0].set_xlim(xmin=0)                                                    # set xmin at 0
                ax1[i,0].set_ylabel('${}$ {}'.format(eigenmotion.columns[i+1], units[i+1]))  # set label of y-axis
                ax1[i,0].minorticks_on()                                                     # set minor ticks
                ax1[i,0].grid(which='major', linestyle='-', linewidth='0.5', color='black')  # customise major grid
                ax1[i,0].grid(which='minor', linestyle=':', linewidth='0.5', color='grey')   # customise minor grid
                ax1[i,0].legend(loc=0, framealpha=1.0).get_frame().set_edgecolor('k')        # set legend for subplot

            ax1[3,0].plot(t, u[0], c='k', linestyle='--', label='Aileron Deflection')        # plot input variable delta_a
            ax1[3,0].plot(t, u[1], c='k', linestyle='-',label='Rudder Deflection')           # plot input variable delta_r
            ax1[3,0].set_xlabel('$t$ [s]')                                                   # set label of y-axis
            ax1[3,0].set_xlim(xmin=0)                                                        # set xmin at 0
            ax1[3,0].set_ylabel('$\delta_a, \delta_r$ [rad]')                                # set label of y-axis
            ax1[3,0].minorticks_on()                                                         # set minor ticks
            ax1[3,0].grid(which='major', linestyle='-', linewidth='0.5', color='black')      # customise major grid
            ax1[3,0].grid(which='minor', linestyle=':', linewidth='0.5', color='grey')       # customise minor grid
            ax1[3,0].legend(loc=0, framealpha=1.0).get_frame().set_edgecolor('k')            # set legend

            fig1.tight_layout(pad=1.0)                                                       # increase spacing between subplots
            # fig1.suptitle('Spiral')                                                          # set title of figure
            fig1.savefig('images/spiralroll.png', dpi=300, bbox_inches='tight')              # save figure

# # plt.show()
