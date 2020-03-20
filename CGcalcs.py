#from CGcalc import *
import pandas as pd
from main import data
import matplotlib.pyplot as plt

time = data['time']
ff_le = data['rh_engine_FMF']
ff_re = data['lh_engine_FMF']
fuel_moment = pd.read_excel('FuelCG.xlsx', header=None, sheet_name='Sheet1')

#-------------functions-------------------------

def CGshift1(Mi, CGi, Mf, mom_change):
    return (Mi*CGi+ mom_change)/(Mf)


def CGshift2(M, CGi, mom_change):
    return CGi + mom_change/M

def CG_MAC(CG_datum,LEMAC, MAC ):
    CG = (CG_datum - LEMAC)/MAC *100
    return CG

def interpolatefuel(fuel, fuel_moment):
      
    for i in range(len(fuel_moment)):
        if fuel<= fuel_moment.iat[i,0]:
            momentcg = fuel_moment.iat[i,1] - ((fuel_moment.iat[i,1] - fuel_moment.iat[i-1,1])/100 )* (fuel_moment.iat[i,0] - fuel)
            break
    return momentcg

def TrapArea(j, time, ff_le, ff_re):
    dt = time[j+1] -time[j]
    f1 = ff_le[j]/3600
    f2 = ff_le[j+1]/3600
    f3 = ff_re[j]/3600
    f4 = ff_re[j+1]/3600
    I1 = (f1+f2)*dt/2
    I2 = (f3+f4)*dt/2
    return I1+I2

def fuelUsed(time, ff_le, ff_re):
    fuelUsed = {}
    integral = 0
    for t in time[:len(time)-1]:
        integral =integral + TrapArea(time[time == t].index[0], time, ff_le, ff_re)
        fuelUsed[t] = integral
    return fuelUsed

def CG_time(t, fuel_i, ZFM, CG_ZFM, fuelUsed):
    fuel_mass = fuel_i - fuelUsed[t]
    fuel_mom = interpolatefuel(fuel_mass, fuel_moment)
    tot_mass = ZFM+ fuel_mass
    CG = CGshift1(ZFM, CG_ZFM, tot_mass, fuel_mom)
    
    return CG, tot_mass
#-------------------------------
MAC = 2.0569 #[m]
LEMAC = 6.64083 #[m]


#--Basic empty mass (BEM) taken from weight measurements
BEM = 9165 #[lbs]
CG_BEM = 291.647954 #inches


#--Payload weight using passenger/payload data
x_seats = [131,131,214,214,251,251,288,288,170]
M_seats = [90*2.20462,102*2.20462,83*2.20462,94*2.20462,84*2.20462,74*2.20462,79*2.20462,103*2.20462,80*2.20462]

x_bag = [74, 321, 338]
M_bag = [0, 0, 0]

M_PL = sum(M_seats) + sum(M_bag)

mom_PL  = sum([x_seats[i]*M_seats[i] for i in range(len(x_seats))]) + sum([x_bag[j]*M_bag[j] for j in range(len(x_bag))])

#ZFW CG:
ZFM = BEM+M_PL
CG_ZFM = CGshift1(BEM, CG_BEM, ZFM, mom_PL) 

#calculate ramp weight CG:
#--Using interpolated fuel-moment chart

Initial_fuel = 4100 #[lbs]

mom_fuel = interpolatefuel(Initial_fuel, fuel_moment)

RM = ZFM + Initial_fuel

CG_RM = CGshift1(ZFM, CG_ZFM, RM, mom_fuel)



#---Plotting----

if __name__ == '__main__':
    
    cgg = []
    m = []
    fuel_used =  fuelUsed(time, ff_le, ff_re)
    for t in time[:len(time)-1]:
        cg, mass = CG_time(t, Initial_fuel, ZFM,CG_ZFM,fuel_used)
        cgg.append(cg)
        m.append(mass)
    plt.figure()
    plt.plot(time[:len(time)-1], m)
    
    plt.figure()
    plt.plot(time[:len(time)-1],cgg)
    
    plt.show()








