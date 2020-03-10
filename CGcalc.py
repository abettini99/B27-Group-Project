import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from main import data

df = pd.read_excel('FuelCG.xlsx', header=None, sheet_name='Sheet1')
time_flight = data['time']
ff_le = data['lh_engine_FMF']
ff_re = data['rh_engine_FMF']
def TrapArea(j):
    dt = time_flight[j+1] -time_flight[j]
    f1 = ff_le[j]/3600
    f2 = ff_le[j+1]/3600
    f3 = ff_re[j]/3600
    f4 = ff_re[j+1]/3600
    I1 = (f2-f1)* (dt)/2. +f1*(dt)
    I2 = (f4-f3)* (dt)/2. +f3*(dt)
    return I1+I2

def totalfuelused(time):
    j = 0    
    
    totalfuelused = 0
    while  time_flight[j]< time :
        
        integral_dt = TrapArea(j)
        totalfuelused = totalfuelused + integral_dt
        j = j + 1
        integral_dt = 0
    
    return totalfuelused


fuelUsed = {}
integral = 0
for t in time_flight[0:48320]:
    integral += TrapArea(time_flight[time_flight == t].index[0])
    fuelUsed[t] = integral


def interpolatefuel(fuel):
    
    for i in range(len(df)):
        if fuel<= df.iat[i,0]:
            momentcg = df.iat[i,1] - ((df.iat[i,1] - df.iat[i-1,1])/100 )* (df.iat[i,0] - fuel)
            break
    
    
    return momentcg


