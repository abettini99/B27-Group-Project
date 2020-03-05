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
