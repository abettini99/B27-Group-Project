import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import math as m
from main import *

TAS = data.Dadc1_tas
u_cap = (TAS[1:] - TAS[:-1])

# speed_of_sound = np.sqrt(1.4 * 287.05 * (data.Dadc1_sat + 273.15))

TAS = manouvre("phugoid").Dadc1_tas

def validation_plot(data_type, flightmanouvre):
    
    if data_type == "u_cap":
        VT0 = np.average(manouvre(flightmanouvre).Dadc1_tas)
        plot_data = (manouvre(flightmanouvre).Dadc1_tas - VT0) / VT0     
    elif data_type == "alpha":
        plot_data = np.radians(manouvre(flightmanouvre).vane_TAS - manouvre(flightmanouvre).vane_TAS[0])     
    elif data_type == "theta":
        plot_data = np.radians(manouvre(flightmanouvre).Ahrs1_Roll - manouvre(flightmanouvre).Ahrs1_Roll[0])   
    elif data_type == "beta":
        



    return plot_data