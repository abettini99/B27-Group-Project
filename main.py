"""
Institution:    TU Delft
Authors:        B27
Date:           06-03-2020

Dynamic Analysis of Citation II
"""

# ==============================================================================================
# Import Libraries
# ==============================================================================================
import pandas as pd                 # package for improved data analysis through DataFrames, etc...
import numpy as np                  # fundamental package for scientific computing
import control.matlab as ml         # import module to emulate functionality of MATLAB
import control as ctl               # import package for analysis and design of feedback control systems
import matplotlib.pyplot as plt     # package to create visualisations
from scipy.io import loadmat        # loadmat imports a .mat file
from numpy.linalg import inv, eig   # inv computes the inverse of a matrix; eig computes eigenvalues of matrix

# ==============================================================================================
# Import Parameters
# ==============================================================================================
from parameters import *        # import parameters for Citation II from parameters.py

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

def manouvre(flightmanouvre):
    """
        This function slices the dataframe into a smaller dataframe for each flight manouvre with the corresponding start and stop time
        :filename: name of flightmanouvre (phugoid, shortperiodoscillation, heavilydampedmotion, spiral or dutchroll)
        :return: sliced dataframe with each variable in one column
    """
    global data                         # declare imported .mat-data in dataframe format as global variable
    if flightmanouvre == "clcd":
        time_start  = 900
        time_stop   = 1600
        data        = data[(data['time'] >= time_start) & (data['time'] <= time_stop)]
        return data
    if flightmanouvre == "elevatortrim":
        time_start  = 1700
        time_stop   = 2100
        data        = data[(data['time'] >= time_start) & (data['time'] <= time_stop)]
        return data
    if flightmanouvre == "cgshift":
        time_start  = 2200
        time_stop   = 2300
        data        = data[(data['time'] >= time_start) & (data['time'] <= time_stop)]
        return data
    if flightmanouvre == "phugoid":
        time_start  = 2600
        time_stop   = 2760
        data        = data[(data['time'] >= time_start) & (data['time'] <= time_stop)]
        return data
    if flightmanouvre == "shortperiodoscillation":
        time_start  = 2700
        time_stop   = 2880
        data        = data[(data['time'] >= time_start) & (data['time'] <= time_stop)]
        return data
    if flightmanouvre == "dutchroll":
        time_start  = 2800
        time_stop   = 3000
        data        = data[(data['time'] >= time_start) & (data['time'] <= time_stop)]
        return data
    if flightmanouvre == "dutchrollYD":
        time_start  = 2900
        time_stop   = 3060
        data        = data[(data['time'] >= time_start) & (data['time'] <= time_stop)]
        return data
    if flightmanouvre == "aperroll":
        time_start  = 3000
        time_stop   = 3240
        data        = data[(data['time'] >= time_start) & (data['time'] <= time_stop)]
        return data
    if flightmanouvre == "spiral":
        time_start  = 3200
        time_stop   = 3480
        data        = data[(data['time'] >= time_start) & (data['time'] <= time_stop)]
        return data

# ==============================================================================================
# Import flight test data from .mat-file
# ==============================================================================================
data    = importdata('referencedata.mat')

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
V       = 1     # [m/s] magnitude of airspeed vector

As[0,0] = - 2 * muc * c / V

As[1,1] = (CZadot - 2 * muc) * c / V

As[2,2] = -V / c

As[3,1] = Cmadot * c / V
As[3,3] = 2 * muc * KY2 * c / V

Bs[0,0] = -CXu
Bs[0,1] = -CXa
Bs[0,2] = -CZ0
Bs[0,3] = -CXq

Bs[1,0] = -CZu
Bs[1,1] = -CZa
Bs[1,2] = -CX0
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
sys = ml.ss(A, B, C, D)         # create state-space system
dt  = np.arange(0, 1000, 0.1)   # create time vector with 0.1s step size

# ==============================================================================================
# Eigenvalue Analysis of matrix A
# ==============================================================================================
ev = eig(A)         # compute eigenvalues and eigenvectors of square matrix A
evals = ev[0]       # eigenvalues of matrix A
evecs = ev[1]       # eigenvectors of matrix A

print('=================== EIGENVALUES OF MATRIX A ===================')
print(evals)
print('===============================================================')

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
i = 0                                                           # running variable for different inputs [de, dt, da, dr]
for df in (initial_de, initial_dt, initial_da, initial_dr):     # iterate over all four lists
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
    # t, y, x = ctl.forced_response(sys, dt, X0)                  # calculate forced response
    # df2 = pd.DataFrame(np.transpose(y), columns=columns)        # convert forced response to DataFrame
    # df.append(df2)                                              # append DataFrame to individual list

# concatenate list into panda dataframe along axis 1
# forced_de = pd.concat(forced_de, axis=1)
# forced_dt = pd.concat(forced_dt, axis=1)
# forced_da = pd.concat(forced_da, axis=1)
# forced_dr = pd.concat(forced_dr, axis=1)
# ==============================================================================================
# FORCED RESPONSE MUST BE FIXED WITH CORRECT INITIAL CONDITION ARRAY X0
# ==============================================================================================

# ==============================================================================================
# Plot step, impulse, initial response of state-space system
# ==============================================================================================
input1    = r'\delta_e'
fig1, ax1 = plt.subplots(2,2, squeeze=False, figsize=(16,9))                                # initialise figure 4 with a (2 x 2) plot layout
for df in (step_de, impulse_de): #, initial_de, initial_dt):
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

# input2 = r'\delta_t'
# fig2, ax2 = plt.subplots(2,2, squeeze=False, figsize=(16,9))                                # initialise figure 4 with a (2 x 2) plot layout
# for df in (step_dt, impulse_dt): #, initial_dt):
#     df = df.loc[:, (df != 0).any(axis=0)]                                                   # remove zero columns for automated plotted
#     ax2[0,0].plot(t, df.iloc[:,0], label='${}$ for ${}$'.format(df.columns[0], input2))     # plot first column in top left plot
#     ax2[0,1].plot(t, df.iloc[:,1], label='${}$ for ${}$'.format(df.columns[1], input2))     # plot second column in top right plot
#     ax2[1,0].plot(t, df.iloc[:,2], label='${}$ for ${}$'.format(df.columns[2], input2))     # plot third column in bottom left plot
#     ax2[1,1].plot(t, df.iloc[:,3], label='${}$ for ${}$'.format(df.columns[3], input2))     # plot fourth column in bottom right plot

#     # Add legends to each subplot
#     ax2[0,0].legend(loc=0, framealpha=1.0).get_frame().set_edgecolor('k')                   # set legend for subplot (0,0)
#     ax2[0,0].set_xlabel('$t$ [s]')                                                          # set label of x-axis for subplot (0,0)
#     ax2[0,0].set_ylabel('$y$ [-]')                                                          # set label of y-axis for subplot (0,0)
#     ax2[0,1].legend(loc=0, framealpha=1.0).get_frame().set_edgecolor('k')                   # set legend for subplot (0,1)
#     ax2[0,1].set_xlabel('$t$ [s]')                                                          # set label of x-axis for subplot (0,1)
#     ax2[0,1].set_ylabel('$y$ [deg]')                                                        # set label of y-axis for subplot (0,1)
#     ax2[1,0].legend(loc=0, framealpha=1.0).get_frame().set_edgecolor('k')                   # set legend for subplot (1,0)
#     ax2[1,0].set_xlabel('$t$ [s]')                                                          # set label of x-axis for subplot (1,0)
#     ax2[1,0].set_ylabel('$y$ [deg]')                                                        # set label of y-axis for subplot (1,0)
#     ax2[1,1].legend(loc=0, framealpha=1.0).get_frame().set_edgecolor('k')                   # set legend for subplot (1,1)
#     ax2[1,1].set_xlabel('$t$ [s]')                                                          # set label of x-axis for subplot (1,1)
#     ax2[1,1].set_ylabel('$y$ [-]')                                                          # set label of y-axis for subplot (1,1)

input3 = r'\delta_a'
fig3, ax3 = plt.subplots(2,2,squeeze=False,figsize=(16,9))                                  # initialise figure 3 with a (2 x 2) plot layout
for df in (step_da, impulse_da): #, initial_da):
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

input4 = r'\delta_r'
fig4, ax4 = plt.subplots(2,2,squeeze=False,figsize=(16,9))                                  # initialise figure 4 with a (2 x 2) plot layout
for df in (step_da, impulse_da): #, initial_da):
    df = df.loc[:, (df != 0).any(axis=0)]                                                   # remove zero columns for automated plotted
    ax4[0,0].plot(t, df.iloc[:,0], label='${}$ for ${}$'.format(df.columns[0], input4))     # plot first column in top left plot
    ax4[0,1].plot(t, df.iloc[:,1], label='${}$ for ${}$'.format(df.columns[1], input4))     # plot second column in top right plot
    ax4[1,0].plot(t, df.iloc[:,2], label='${}$ for ${}$'.format(df.columns[2], input4))     # plot third column in bottom left plot
    ax4[1,1].plot(t, df.iloc[:,3], label='${}$ for ${}$'.format(df.columns[3], input4))     # plot fourth column in bottom right plot

    # Add legends to each subplot
    ax4[0,0].legend(loc=0, framealpha=1.0).get_frame().set_edgecolor('k')                   # set legend for subplot (0,0)
    ax4[0,0].set_xlabel('$t$ [s]')                                                          # set label of x-axis for subplot (0,0)
    ax4[0,0].set_ylabel('$y$ [-]')                                                          # set label of y-axis for subplot (0,0)
    ax4[0,1].legend(loc=0, framealpha=1.0).get_frame().set_edgecolor('k')                   # set legend for subplot (0,1)
    ax4[0,1].set_xlabel('$t$ [s]')                                                          # set label of x-axis for subplot (0,1)
    ax4[0,1].set_ylabel('$y$ [deg]')                                                        # set label of y-axis for subplot (0,1)
    ax4[1,0].legend(loc=0, framealpha=1.0).get_frame().set_edgecolor('k')                   # set legend for subplot (1,0)
    ax4[1,0].set_xlabel('$t$ [s]')                                                          # set label of x-axis for subplot (1,0)
    ax4[1,0].set_ylabel('$y$ [deg]')                                                        # set label of y-axis for subplot (1,0)
    ax4[1,1].legend(loc=0, framealpha=1.0).get_frame().set_edgecolor('k')                   # set legend for subplot (1,1)
    ax4[1,1].set_xlabel('$t$ [s]')                                                          # set label of x-axis for subplot (1,1)
    ax4[1,1].set_ylabel('$y$ [-]')                                                          # set label of y-axis for subplot (1,1)

# Save figures for each input variable
fig1.savefig('images/response_de.png', dpi=300, bbox_inches='tight')
# fig2.savefig('images/response_dt.png', dpi=300, bbox_inches='tight')
fig3.savefig('images/response_da.png', dpi=300, bbox_inches='tight')
fig4.savefig('images/response_dr.png', dpi=300, bbox_inches='tight')

plt.show()





