"""
Institution:    TU Delft
Authors:        B27
Date:           06-03-2020

Dynamic Analysis of Citation II
"""

## =========== Library Imports          ===========
from scipy.io import loadmat    # loadmat imports a .mat file
from numpy.linalg import inv    # inv computes the inverse of a matrix
import pandas as pd
import numpy as np

## =========== Parameter Imports        ===========
from parameters import *                # import parameters for Citation II from parameters.py

## =========== Function Definitions     ===========
def importdata(filename):
    """
        This function imports the data from the flight tests
        :filename: relative path of the .mat-file from the flight test
        :return: dataframe with each variable in one column
    """

    mat = loadmat(filename) # load data from .mat-file - returns dictionary with variable names as keys and loaded matrices as values
    mdata = mat['flightdata'] # access first level variable in .mat-file
    mdtype = mdata.dtype # dtypes of structures are unsized objects
    ndata = {n: mdata[n][0,0]['data'][0,0] for n in mdtype.names} # redefine data as dictionary with names from dtypes

    data = pd.DataFrame(data=None, index=None, columns=None, dtype=None, copy=False) # initialise empty dataframe

    for key, values in ndata.items(): # iterate over keys and
        data[key] = values.flatten() # add new column with variable as key and values; input to dataframe must be 1D such that 2D arrays must be flattened

    return data

def manouvre(flightmanouvre):
    """
        This function slices the dataframe into a smaller dataframe for each flight manouvre with the corresponding start and stop time
        :filename: name of flightmanouvre (phugoid, shortperiodoscillation, heavilydampedmotion, spiral or dutchroll)
        :return: sliced dataframe with each variable in one column
    """

    global data # declare imported .mat-data in dataframe format as global variable
    if flightmanouvre == "phugoid":
        time_start = 20
        time_stop = 50
        data = data[(data['time'] >= time_start) & (data['time'] <= time_stop)]
        return data
    if flightmanouvre == "shortperiodoscillation":
        time_start = 20
        time_stop = 50
        data = data[(data['time'] >= time_start) & (data['time'] <= time_stop)]
        return data
    if flightmanouvre == "heavilydampedmotion":
        time_start = 20
        time_stop = 50
        data = data[(data['time'] >= time_start) & (data['time'] <= time_stop)]
        return data
    if flightmanouvre == "spiral":
        time_start = 20
        time_stop = 50
        data = data[(data['time'] >= time_start) & (data['time'] <= time_stop)]
        return data
    if flightmanouvre == "dutchroll":
        time_start = 20
        time_stop = 50
        data = data[(data['time'] >= time_start) & (data['time'] <= time_stop)]
        return data

## ============== Main File ==============
# Import flight test data from .mat-file
data = importdata('referencedata.mat')
#print("\n")
# print(data)


# Declaration of matrices and column vectors
A = np.zeros((8,8))         # Declaration of matrix A with dimensions [8 x 8] for system of equations
B = np.zeros((8,8))         # Declaration of matrix B with dimensions [8 x 4] for system of equations
C = np.zeros((8,4))         # Declaration of matrix C with dimensions [8 x 4] for system of equations
A_temp = np.zeros((8,8))    # Declaration of temporary matrix A with dimensions [8 x 8] for system of equations
B_temp = np.zeros((8,8))    # Declaration of temporary matrix B with dimensions [8 x 4] for system of equations
C_temp = np.zeros((8,4))    # Declaration of temporary matrix C with dimensions [8 x 4] for system of equations
As = np.zeros((4,4))        # Declaration of matrix As with dimensions [4 x 4] for symmetric EOM
Aa = np.zeros((4,4))        # Declaration of matrix Aa with dimensions [4 x 4] for asymmetric EOM
Bs = np.zeros((4,4))        # Declaration of matrix Bs with dimensions [4 x 2] for symmetric EOM
Ba = np.zeros((4,4))        # Declaration of matrix Ba with dimensions [4 x 2] for asymmetric EOM
Cs = np.zeros((4,2))        # Declaration of matrix Cs with dimensions [4 x 2] for symmetric EOM
Ca = np.zeros((4,2))        # Declaration of matrix Ca with dimensions [4 x 2] for asymmetric EOM

# FOLLOWING VARIABLES ARE NOT DEFINED IN PARAMETERS.PY
V = 1
CXdt = 1
CZdt = 1
Cmdt = 1

# Population of symmetric EOM matrices with variables for state-space representation
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
Cs[0,1] = -CXdt

Cs[1,0] = -CZde
Cs[1,1] = -CZdt

Cs[3,0] = -Cmde
Cs[3,1] = -Cmdt

# Population of asymmetric EOM matrices with variables for state-space representation
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

# Population of matrix A with matrices As and Aa
A_temp[0:4, 0:4] = As
A_temp[4:8, 4:8] = Aa
B_temp[0:4, 0:4] = Bs
B_temp[4:8, 4:8] = Ba
C_temp[0:4, 0:2] = Cs
C_temp[4:8, 2:4] = Ca

# Compute final matrices A and B based equation given in simulation plan
A = np.dot(inv(A_temp), B_temp)
B = np.dot(inv(A_temp), C_temp)
