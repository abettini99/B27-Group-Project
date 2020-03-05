from scipy.io import loadmat
import pandas as pd
import numpy as np

# import parameters for Citation II from parameters.py
from parameters import *


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

# Import flight test data from .mat-file
data = importdata('referencedata.mat')


# Declaration of matrices and column vectors
A = np.zeros((8,8))     # Declaration of matrix A with dimensions [8 x 8] for system of equations
B = np.zeros((8,4))     # Declaration of matrix B with dimensions [8 x 4] for system of equations
As = np.zeros((4,4))    # Declaration of matrix As with dimensions [4 x 4] for symmetric EOM
Aa = np.zeros((4,4))    # Declaration of matrix Aa with dimensions [4 x 4] for asymmetric EOM
Bs = np.zeros((4,2))    # Declaration of matrix Bs with dimensions [4 x 2] for symmetric EOM
Ba = np.zeros((4,2))    # Declaration of matrix Ba with dimensions [4 x 2] for asymmetric EOM

# FOLLOWING VARIABLES ARE NOT DEFINED IN PARAMETERS.PY
V = 1
CXdt = 1
CZdt = 1
Cmdt = 1

# Population of symmetric EOM matrices with variables for state-space representation
k1  = V / c * 1 / (2 * muc)
k2  = V / c * 1 / (2 * muc - CZadot)
k3  = V / c * 1 / (2 * muc * KY2)
k4  = 1 / (2 * muc - CZadot)

As[0,0] = k1 * CXu
As[0,1] = k1 * CXa
As[0,2] = k1 * CZ0

As[1,0] = k2 * CZu
As[1,1] = k2 * CZa
As[1,2] = -k2 * CX0
As[1,3] = k2 * (2 * muc + CZq)

As[2,3] = V / c

As[3,0] = k3  * (Cmu + CZu * Cmadot * k4)
As[3,1] = k3 * (Cma + CZa * Cmadot * k4)
As[3,2] = -k3 * CX0 * Cmadot * k4
As[3,3] = k3 * (Cmq + Cmadot * (2 * muc + CZq) / k4)

Bs[0,0] = k1 * CXde
Bs[0,1] = k1 * CXdt

Bs[1,0] = k2 * CZde
Bs[1,1] = k2 * CZdt

Bs[3,0] = k3 * (Cmde + CZde * Cmadot * k4)
Bs[3,1] = k3 * (Cmdt + CZdt * Cmadot * k4)

# Population of asymmetric EOM matrices with variables for state-space representation
k5 = V / b * 1 / (2 * mub)
k6 = V / b * 1 / (4 * mub * (KX2 * KZ2 - KXZ))

Aa[0,0] = k5 * CYb
Aa[0,1] = k5 * CL
Aa[0,2] = k5 * CYp
Aa[0,3] = k5 * (CYr - 4 * mub)

Aa[1,2] = 2 * V / b

Aa[2,0] = k6 * (Clb * KZ2 + Cnb * KXZ)
Aa[2,2] = k6 * (Clp * KZ2 + Cnp * KXZ)
Aa[2,3] = k6 * (Clr * KZ2 + Cnr * KXZ)

Aa[3,0] = k6 * (Clb * KXZ + Cnb * KX2)
Aa[3,2] = k6 * (Clp * KXZ + Cnp * KX2)
Aa[3,3] = k6 * (Clr * KXZ + Cnr * KX2)

Ba[0,1] = k5 * CYdr

Ba[2,0] = k6 * (Clda * KZ2 + Cnda * KXZ)
Ba[2,1] = k6 * (Cldr * KZ2 + Cndr * KXZ)

Ba[3,0] = k6 * (Clda * KXZ + Cnda * KX2)
Ba[3,1] = k6 * (Cldr * KXZ + Cndr * KX2)

# Population of matrix A with matrices As and Aa
A[0:4, 0:4] = As
A[4:8, 4:8] = Aa
B[0:4, 0:2] = Bs
B[4:8, 2:8] = Ba
