import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


df = pd.read_excel('FuelCG.xlsx', header=None, sheet_name='Sheet1')

def integralfuel(j):
    dt = 0.1
    f1 = 3
    f2 = 1
    I = (f2-f1)* (dt)/2. +f1*(dt)
    return I

def integral_time(time):
    j = 0    
    dt = 0.1
    integral_dt = 0
    while  j< time :
        
        integral_dt = integralfuel(j)
        integral_time = integral_time + integral_dt
        j = j + dt
        integralfuel = 0
    
    return integral_time


def interpolatefuel(fuel):
    
    for i in range(len(df)):
        if fuel<= df.iat[i,0]:
            momentcg = df.iat[i,1] - ((df.iat[i,1] - df.iat[i-1,1])/100 )* (df.iat[i,0] - fuel)
            break
    
    cg  = momentcg*100/fuel
    return cg


