"""
Institution:    TU Delft
Authors:        B27
Date:           06-03-2020

Dynamic Analysis of Citation II
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
        time_stop   = 2820
        data        = data[(data['time'] >= time_start) & (data['time'] <= time_stop)]
        return data

    if flightmanouvre == "shortperiod":
        time_start  = 2820 - 10
        time_stop   = 2880
        data        = data[(data['time'] >= time_start) & (data['time'] <= time_stop)]
        return data

    if flightmanouvre == "dutchroll":
        time_start  = 2880 - 10
        time_stop   = 3000
        data        = data[(data['time'] >= time_start) & (data['time'] <= time_stop)]
        return data
    if flightmanouvre == "dutchrollYD":
        time_start  = 3000 - 10
        time_stop   = 3060
        data        = data[(data['time'] >= time_start) & (data['time'] <= time_stop)]
        return data

    if flightmanouvre == "aperroll":
        time_start  = 3060 - 10
        time_stop   = 3240
        data        = data[(data['time'] >= time_start) & (data['time'] <= time_stop)]
        return data

    if flightmanouvre == "spiral":
        time_start  = 3240 - 10
        time_stop   = 3480
        data        = data[(data['time'] >= time_start) & (data['time'] <= time_stop)]
        return data

def ktstoms(velocity):
    """
        Function converts velocity given in knots [kts] to metre per second [m/s]
        :velocity: input velocity in knots
        :return: velocity converted to metre per second
    """
    return velocity * 0.5144444444

def fttom(altitude):
    """
        Function converts altitude given in feet [ft] to metre [m]
        :velocity: input altitude in [ft]
        :return: outputs altitude in [m]
    """
    return altitude * 0.3048

# ==============================================================================================
# Import data from Matlab files and transform coordinate system from body to stability axis
# ==============================================================================================
# data = importdata('referencedata.mat')  # initialise reference data from matlab file
data = importdata('flightdata.mat')     # initialise flight data from matlab file

alpha0 = radians(data.vane_AOA.iloc[9830])    # [rad] angle of attack in the stationary flight condition
theta0 = radians(data.Ahrs1_Pitch.iloc[9830]) # [rad] pitch angle in the stationary flight condition

# Transform angle of attack and pitch angle from body to stability axis frame
data['vane_AOA'] = data['vane_AOA'] - degrees(alpha0)
data['Ahrs1_Pitch'] = data['Ahrs1_Pitch'] - degrees(theta0)

# ==============================================================================================
# Stationary measurements
# ==============================================================================================
# data    = importdata('referencedata.mat')   # initialise flight data from matlab file
# clcd    = manouvre('clcd')                  # sliced data for the 6 CL-CD measurement series
# etrim   = manouvre('elevatortrim')          # sliced data for the 7 e-trim measurement series
# cgshift = manouvre('cgshift')               # sliced data for the 2 cg-shift measurement series

# ==============================================================================================
# Eigenmotion analysis - uncomment required eigenmotion array
# ==============================================================================================
# data = manouvre(data, 'phugoid')                     # sliced data array for phugoid motion
# data = manouvre(data, 'shortperiod')                 # sliced data array short period oscillation motion
data = manouvre(data, 'dutchroll')                   # sliced data array for dutch roll motion
# data = manouvre(data, 'dutchrollYD')                 # sliced data array for yawed dutch roll motion
# data = manouvre(data, 'aperroll')                    # sliced data array for aperiodic roll motion
# data = manouvre(data, 'spiral')                      # sliced data array for spiral motion

# ==============================================================================================
# Parameter definition; copied from Cit_par.py
# ==============================================================================================
hp0    = fttom(data.Dadc1_alt.iloc[0:10].mean())     # [m] pressure altitude in the stationary flight condition
V0     = ktstoms(data.Dadc1_tas.iloc[0:10].mean())   # [m/s] true airspeed in the stationary flight condition
alpha0 = radians(data.vane_AOA.iloc[0:10].mean())    # [rad] angle of attack in the stationary flight condition
theta0 = radians(data.Ahrs1_Pitch.iloc[0:10].mean()) # [rad] pitch angle in the stationary flight condition

m      = 6805.903           # [kg] takeoff weight of Cessna Citation II

e      = 0.8                # [-] Oswald factor
CD0    = 0.04               # [-] Zero lift drag coefficient
CLa    = 5.084              # [-] Slope of CL-alpha curve

# MUST BE FILLED IN BASED ON SECOND MEASUREMENT SERIES - CURRENTLY SET FURTHER DOWN FROM REFERENCE VALUES
Cma    = -0.582128          # [-] longitudinal stabilty
Cmde   = -1.21076           # [-] elevator effectiveness - second equation after (7-26) in FD reader

S      = 30.00              # [m^2] wing area
Sh     = 0.2 * S            # [m^2] stabiliser area
Sh_S   = Sh / S             # [-]
lh     = 0.71 * 5.968       # [m] tail length
c      = 2.0569             # [m] mean aerodynamic cord
lh_c   = lh / c             # [-]
b      = 15.911             # [m] wing span
bh     = 5.791              # [m] stabilser span
AR     = b ** 2 / S         # [-] wing aspect ratio
Ah     = bh ** 2 / Sh       # [-] stabilser aspect ratio
Vh_V   = 1                  # [-]
ih     = -2 * pi / 180      # [rad] stabiliser angle of incidence

rho0   = 1.2250             # [kg/m^3] air density at sea level
LAMBDA = -0.0065            # [K/m] temperature gradient in ISA
Temp0  = 288.15             # [K] temperature at sea level in ISA
R      = 287.05             # [m^2/s^2 K] specific gas constant
g      = 9.81               # [m/s^2] gravitational acceleration
rho    = rho0 * pow( ((1+(LAMBDA * hp0 / Temp0))), (-((g / (LAMBDA*R)) + 1))) # [kg/m^3] density at altitude h
W      = m * g              # [N] aircraft weight

muc    = m / (rho * S * c)
mub    = m / (rho * S * b)
KX2    = 0.019
KY2    = 1.3925
KZ2    = 0.042
KXZ    = 0.002

Cmac   = 0                                      # [-] Moment coefficient about the aerodynamic centre
CNwa   = CLa                                    # [-] Wing normal force slope
CNha   = 2 * pi * Ah / (Ah + 2)                 # [-] Stabiliser normal force slope
depsda = 4 / (AR + 2)                           # [-] Downw sh gradient

CL     = 2 * W / (rho * V0 ** 2 * S)                # [-] Lift coefficient
CD     = CD0 + (CLa * alpha0) ** 2 / (pi * AR * e)  # [-] Drag coefficient

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
CZde   = -0.6961

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
A       = np.zeros((8,8))         # Declaration of matrix A with dimensions [8 x 8] for system of equations
B       = np.zeros((8,8))         # Declaration of matrix B with dimensions [8 x 4] for system of equations
C       = np.zeros((8,4))         # Declaration of matrix C with dimensions [8 x 4] for system of equations
A_temp  = np.zeros((8,8))         # Declaration of temporary matrix A with dimensions [8 x 8] for system of equations
B_temp  = np.zeros((8,8))         # Declaration of temporary matrix B with dimensions [8 x 4] for system of equations
C_temp  = np.zeros((8,4))         # Declaration of temporary matrix C with dimensions [8 x 4] for system of equations
As      = np.zeros((4,4))         # Declaration of matrix As with dimensions [4 x 4] for symmetric EOM
Aa      = np.zeros((4,4))         # Declaration of matrix Aa with dimensions [4 x 4] for asymmetric EOM
Bs      = np.zeros((4,4))         # Declaration of matrix Bs with dimensions [4 x 2] for symmetric EOM
Ba      = np.zeros((4,4))         # Declaration of matrix Ba with dimensions [4 x 2] for asymmetric EOM
Cs      = np.zeros((4,2))         # Declaration of matrix Cs with dimensions [4 x 2] for symmetric EOM
Ca      = np.zeros((4,2))         # Declaration of matrix Ca with dimensions [4 x 2] for asymmetric EOM

# ==============================================================================================
# Population of symmetric EOM matrices with variables for state-space representation
# ==============================================================================================
V       = V0    # [m/s] redefine magnitude of airspeed vector as true airspeed

As[0,0] = - 2 * muc * c / V

As[1,1] = (CZadot - 2 * muc) * c / V

As[2,2] = -c / V

As[3,1] = Cmadot * c / V
As[3,3] = -2 * muc * KY2 * c / V

Bs[0,0] = -CXu
Bs[0,1] = -CXa
Bs[0,2] = -CZ0
Bs[0,3] = -CXq

Bs[1,0] = -CZu
Bs[1,1] = -CZa
Bs[1,2] = CX0
Bs[1,3] = -(CZq + 2 * muc)

Bs[2,3] = -1

Bs[3,0] = -Cmu
Bs[3,1] = -Cma
Bs[3,3] = -Cmq

Cs[0,0] = -CXde

Cs[1,0] = -CZde

Cs[3,0] = -Cmde

# ==============================================================================================
# Population of asymmetric EOM matrices with variables for state-space representation
# ==============================================================================================
Aa[0,0] = (CYbdot - 2 * mub) * b / V

Aa[1,1] = -0.5 * b / V

Aa[2,2] = -4 * mub * KX2 * b / V
Aa[2,3] = 4 * mub * KXZ * b / V

Aa[3,0] = Cnbdot * b / V
Aa[3,2] = 4 * mub * KXZ * b / V
Aa[3,3] = -4 * mub * KZ2 * b / V

Ba[0,0] = -CYb
Ba[0,1] = -CL
Ba[0,2] = -CYp
Ba[0,3] = -(CYr - 4 * mub)

Ba[1,2] = -1

Ba[2,0] = -Clb
Ba[2,2] = -Clp
Ba[2,3] = -Clr

Ba[3,0] = -Cnb
Ba[3,2] = -Cnp
Ba[3,3] = -Cnr

Ca[0,0] = -CYda
Ca[0,1] = -CYdr

Ca[2,0] = -Clda
Ca[2,1] = -Cldr

Ca[3,0] = -Cnda
Ca[3,1] = -Cndr

# ==============================================================================================
# Population of matrix A with matrices As and Aa
# ==============================================================================================
A_temp[0:4, 0:4] = As
A_temp[4:8, 4:8] = Aa
B_temp[0:4, 0:4] = Bs
B_temp[4:8, 4:8] = Ba
C_temp[0:4, 0:2] = Cs
C_temp[4:8, 2:4] = Ca

# ==============================================================================================
# Compute final matrices A and B based equation given in simulation plan
# ==============================================================================================
A = np.dot(inv(A_temp), B_temp)
B = np.dot(inv(A_temp), C_temp)

# Output of state-space representation should be equal to the relevant aircraft states
# --> matrix C is the identity matrix and D is a zero array
C = np.identity(8)
D = np.zeros((8,4))

# Calculate state-space representation of system for different responses
sys = ml.ss(A, B, C, D)       # create state-space system
dt  = np.arange(0, 50, 0.1)   # create time vector with 0.1s step size

# ==============================================================================================
# Eigenvalue Analysis of matrix A
# ==============================================================================================
evals, evecs = eig(A)         # compute eigenvalues and eigenvectors of square matrix A

print('=================== EIGENVALUES OF MATRIX A ===================')
print(evals)
print('===============================================================')

# # Plot eigenvalues
# fig, ax = plt.subplots(1)
# x = [ev.real for ev in evals]
# y = [ev.imag for ev in evals]
# ax.scatter(x, y, c='r')
# ax.set_xlabel('Re($\lambda$)')
# ax.set_ylabel('Im($\lambda$)')
# plt.grid()
# plt.show()

# ==============================================================================================
# Calculates responses to state-space system
# ==============================================================================================
columns = [r'\hat{u}_s', r'\alpha_s', r'\theta_s', r'\frac{qc}{V}_s', r'\beta_s', r'\phi_s', r'\frac{pb}{2V}_s', r'\frac{rb}{2V}_s']     # names of invidiual columns for DataFrame
step_de, step_dt, step_da, step_dr = [], [], [], []             # initialise lists for step reponse for all four inputs
X0  = np.zeros((8,1))                                           # initial condition for step response
i   = 0                                                         # running variable for different inputs [de, dt, da, dr]
for df in (step_de, step_dt, step_da, step_dr):                 # iterate over all four lists
    t, y  = ctl.step_response(sys, dt, X0, input=i)             # calculate step response
    df2   = pd.DataFrame(np.transpose(y), columns=columns)      # convert step response to DataFrame
    df.append(df2)                                              # append DataFrame to individual list
    i += 1                                                      # change integer to change input variable

# concatenate list into panda dataframe along axis 1
step_de = pd.concat(step_de, axis=1)
step_dt = pd.concat(step_dt, axis=1)
step_da = pd.concat(step_da, axis=1)
step_dr = pd.concat(step_dr, axis=1)

columns = [r'\hat{u}_{im}', r'\alpha_{im}', r'\theta_{im}', r'\frac{qc}{V}_{im}', r'\beta_{im}', r'\phi_{im}', r'\frac{pb}{2V}_{im}', r'\frac{rb}{2V}_{im}']     # names of invidiual columns for DataFrame
impulse_de, impulse_dt, impulse_da, impulse_dr = [], [], [], [] # initialise lists for step reponse for all four inputs
i = 0                                                           # running variable for different inputs [de, dt, da, dr]
for df in (impulse_de, impulse_dt, impulse_da, impulse_dr):     # iterate over all four lists
    t, y = ctl.impulse_response(sys, dt, X0, input=i)           # calculate impulse response
    df2  = pd.DataFrame(np.transpose(y), columns=columns)       # convert impulse response to DataFrame
    df.append(df2)                                              # append DataFrame to individual list
    i += 1                                                      # change integer to change input variable

# concatenate list into panda dataframe along axis 1
impulse_de = pd.concat(impulse_de, axis=1)
impulse_dt = pd.concat(impulse_dt, axis=1)
impulse_da = pd.concat(impulse_da, axis=1)
impulse_dr = pd.concat(impulse_dr, axis=1)

columns = [r'\hat{u}_{in}', r'\alpha_{in}', r'\theta_{in}', r'\frac{qc}{V}_{in}', r'\beta_{in}', r'\phi_{in}', r'\frac{pb}{2V}_{in}', r'\frac{rb}{2V}_{in}']     # names of invidiual columns for DataFrame
initial_de, initial_dt, initial_da, initial_dr = [], [], [], [] # initialise lists for step reponse for all four inputs
deflectionlist = [data.delta_e.iloc[0], 0, data.delta_a.iloc[0], data.delta_r.iloc[0]]  # list to iterate over the initial deflection for each response
i = 0                                                           # running variable for different inputs [de, dt, da, dr]
for df in (initial_de, initial_dt, initial_da, initial_dr):     # iterate over all four lists
    X0[:,] = deflectionlist[i]                                  # populate rows of initial condition with initial deflection
    t, y = ctl.initial_response(sys, dt, X0, input=i)           # calculate initial response
    df2  = pd.DataFrame(np.transpose(y), columns=columns)       # convert initial response to DataFrame
    df.append(df2)                                              # append DataFrame to individual list
    i += 1                                                      # change integer to change input variable

# concatenate list into panda dataframe along axis 1
initial_de = pd.concat(initial_de, axis=1)
initial_dt = pd.concat(initial_dt, axis=1)
initial_da = pd.concat(initial_da, axis=1)
initial_dr = pd.concat(initial_dr, axis=1)

# ==============================================================================================
# ==============================================================================================
# FORCED RESPONSE MUST BE FIXED WITH CORRECT INITIAL CONDITION ARRAY X0
# ==============================================================================================
# ==============================================================================================
# forced_de, forced_dt, forced_da, forced_dr = [], [], [], []     # initialise lists for step reponse for all four inputs
# X0 = np.zeros((4,2000))
# for df in (forced_de, forced_dt, forced_da, forced_dr):         # iterate over all four lists
#     t, y, x = ctl.forced_response(sys, dt, X0)                  # calculate forced response
#     df2 = pd.DataFrame(np.transpose(y), columns=columns)        # convert forced response to DataFrame
#     df.append(df2)                                              # append DataFrame to individual list

# concatenate list into panda dataframe along axis 1
# forced_de = pd.concat(forced_de, axis=1)
# forced_dt = pd.concat(forced_dt, axis=1)
# forced_da = pd.concat(forced_da, axis=1)
# forced_dr = pd.concat(forced_dr, axis=1)
# ==============================================================================================
# FORCED RESPONSE MUST BE FIXED WITH CORRECT INITIAL CONDITION ARRAY X0
# ==============================================================================================

# ==============================================================================================
# Plot step, impulse, initial and forced response of state-space system
# ==============================================================================================
input1    = r'\delta_e'
fig1, ax1 = plt.subplots(2,2, squeeze=False, figsize=(16,9))                                # initialise figure 4 with a (2 x 2) plot layout
for df in (step_de, impulse_de, initial_de): #, forced_de):
    df = df.loc[:, (df != 0).any(axis=0)]                                                   # remove zero columns for automated plotted
    ax1[0,0].plot(t, df.iloc[:,0], label='${}$ for ${}$'.format(df.columns[0], input1))     # plot first column in top left plot
    ax1[0,1].plot(t, df.iloc[:,1], label='${}$ for ${}$'.format(df.columns[1], input1))     # plot second column in top right plot
    ax1[1,0].plot(t, df.iloc[:,2], label='${}$ for ${}$'.format(df.columns[2], input1))     # plot third column in bottom left plot
    ax1[1,1].plot(t, df.iloc[:,3], label='${}$ for ${}$'.format(df.columns[3], input1))     # plot fourth column in bottom right plot

    # Add legends to each subplot
    ax1[0,0].legend(loc=0, framealpha=1.0).get_frame().set_edgecolor('k')                   # set legend for subplot (0,0)
    ax1[0,0].set_xlabel('$t$ [s]')                                                          # set label of x-axis for subplot (0,0)
    ax1[0,0].set_ylabel('$y$ [-]')                                                          # set label of y-axis for subplot (0,0)
    ax1[0,1].legend(loc=0, framealpha=1.0).get_frame().set_edgecolor('k')                   # set legend for subplot (0,1)
    ax1[0,1].set_xlabel('$t$ [s]')                                                          # set label of x-axis for subplot (0,1)
    ax1[0,1].set_ylabel('$y$ [deg]')                                                        # set label of y-axis for subplot (0,1)
    ax1[1,0].legend(loc=0, framealpha=1.0).get_frame().set_edgecolor('k')                   # set legend for subplot (1,0)
    ax1[1,0].set_xlabel('$t$ [s]')                                                          # set label of x-axis for subplot (1,0)
    ax1[1,0].set_ylabel('$y$ [deg]')                                                        # set label of y-axis for subplot (1,0)
    ax1[1,1].legend(loc=0, framealpha=1.0).get_frame().set_edgecolor('k')                   # set legend for subplot (1,1)
    ax1[1,1].set_xlabel('$t$ [s]')                                                          # set label of x-axis for subplot (1,1)
    ax1[1,1].set_ylabel('$y$ [-]')                                                          # set label of y-axis for subplot (1,1)

input2    = r'\delta_a'
fig2, ax2 = plt.subplots(2,2,squeeze=False,figsize=(16,9))                                  # initialise figure 3 with a (2 x 2) plot layout
for df in (step_da, impulse_da, initial_da): #, forced_da):
    df = df.loc[:, (df != 0).any(axis=0)]                                                   # remove zero columns for automated plotted
    ax2[0,0].plot(t, df.iloc[:,0], label='${}$ for ${}$'.format(df.columns[0], input2))     # plot first column in top left plot
    ax2[0,1].plot(t, df.iloc[:,1], label='${}$ for ${}$'.format(df.columns[1], input2))     # plot second column in top right plot
    ax2[1,0].plot(t, df.iloc[:,2], label='${}$ for ${}$'.format(df.columns[2], input2))     # plot third column in bottom left plot
    ax2[1,1].plot(t, df.iloc[:,3], label='${}$ for ${}$'.format(df.columns[3], input2))     # plot fourth column in bottom right plot

    # Add legends to each subplot
    ax2[0,0].legend(loc=0, framealpha=1.0).get_frame().set_edgecolor('k')                   # set legend for subplot (0,0)
    ax2[0,0].set_xlabel('$t$ [s]')                                                          # set label of x-axis for subplot (0,0)
    ax2[0,0].set_ylabel('$y$ [-]')                                                          # set label of y-axis for subplot (0,0)
    ax2[0,1].legend(loc=0, framealpha=1.0).get_frame().set_edgecolor('k')                   # set legend for subplot (0,1)
    ax2[0,1].set_xlabel('$t$ [s]')                                                          # set label of x-axis for subplot (0,1)
    ax2[0,1].set_ylabel('$y$ [deg]')                                                        # set label of y-axis for subplot (0,1)
    ax2[1,0].legend(loc=0, framealpha=1.0).get_frame().set_edgecolor('k')                   # set legend for subplot (1,0)
    ax2[1,0].set_xlabel('$t$ [s]')                                                          # set label of x-axis for subplot (1,0)
    ax2[1,0].set_ylabel('$y$ [deg]')                                                        # set label of y-axis for subplot (1,0)
    ax2[1,1].legend(loc=0, framealpha=1.0).get_frame().set_edgecolor('k')                   # set legend for subplot (1,1)
    ax2[1,1].set_xlabel('$t$ [s]')                                                          # set label of x-axis for subplot (1,1)
    ax2[1,1].set_ylabel('$y$ [-]')                                                          # set label of y-axis for subplot (1,1)

input3    = r'\delta_r'
fig3, ax3 = plt.subplots(2,2,squeeze=False,figsize=(16,9))                                  # initialise figure 4 with a (2 x 2) plot layout
for df in (step_dr, impulse_dr, initial_dr): #, forced_dr):
    df = df.loc[:, (df != 0).any(axis=0)]                                                   # remove zero columns for automated plotted
    ax3[0,0].plot(t, df.iloc[:,0], label='${}$ for ${}$'.format(df.columns[0], input3))     # plot first column in top left plot
    ax3[0,1].plot(t, df.iloc[:,1], label='${}$ for ${}$'.format(df.columns[1], input3))     # plot second column in top right plot
    ax3[1,0].plot(t, df.iloc[:,2], label='${}$ for ${}$'.format(df.columns[2], input3))     # plot third column in bottom left plot
    ax3[1,1].plot(t, df.iloc[:,3], label='${}$ for ${}$'.format(df.columns[3], input3))     # plot fourth column in bottom right plot

    # Add legends to each subplot
    ax3[0,0].legend(loc=0, framealpha=1.0).get_frame().set_edgecolor('k')                   # set legend for subplot (0,0)
    ax3[0,0].set_xlabel('$t$ [s]')                                                          # set label of x-axis for subplot (0,0)
    ax3[0,0].set_ylabel('$y$ [-]')                                                          # set label of y-axis for subplot (0,0)
    ax3[0,1].legend(loc=0, framealpha=1.0).get_frame().set_edgecolor('k')                   # set legend for subplot (0,1)
    ax3[0,1].set_xlabel('$t$ [s]')                                                          # set label of x-axis for subplot (0,1)
    ax3[0,1].set_ylabel('$y$ [deg]')                                                        # set label of y-axis for subplot (0,1)
    ax3[1,0].legend(loc=0, framealpha=1.0).get_frame().set_edgecolor('k')                   # set legend for subplot (1,0)
    ax3[1,0].set_xlabel('$t$ [s]')                                                          # set label of x-axis for subplot (1,0)
    ax3[1,0].set_ylabel('$y$ [deg]')                                                        # set label of y-axis for subplot (1,0)
    ax3[1,1].legend(loc=0, framealpha=1.0).get_frame().set_edgecolor('k')                   # set legend for subplot (1,1)
    ax3[1,1].set_xlabel('$t$ [s]')                                                          # set label of x-axis for subplot (1,1)
    ax3[1,1].set_ylabel('$y$ [-]')                                                          # set label of y-axis for subplot (1,1)

# Save figures for each input variable
fig1.savefig('images/response_de.png', dpi=300, bbox_inches='tight')
fig2.savefig('images/response_da.png', dpi=300, bbox_inches='tight')
fig3.savefig('images/response_dr.png', dpi=300, bbox_inches='tight')

# plt.show()
