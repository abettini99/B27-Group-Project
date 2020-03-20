from CGcalcs import *
import numpy as np
import pandas as pd

# UNIT 1
#Add 1000kg at 10m, together with 1000kg initial mass at 5m >>>  final cg of 7.5 is expected.
print(CGshift1(1000, 5, 2000, 10000))

#unit 2
print(CGshift2(1000, 5, -500))


#unit 3:
#160lbs of fuel should have a moment 60% between  298.16 and 591.18  = 473.972
print(interpolatefuel(160, fuel_moment))

#unit 4 (higher level system):
# each engine has fuel flow linearly from 200 to 500, during a time from 0 to 100 with time step of 0.1. Integral is calculated with the trapezoidal rule
#and calculated by hand by calculating the area under the graph. Result is 19.4444. This should be exact with no numerical error.
time2 = np.arange(0,100.1,0.1)
time2 = pd.Series(time2)
ff = np.linspace(200,500,len(time2))

fuelused = fuelUsed(time2,ff,ff)
print(fuelused) 



#System test 
cgg = []
m = []
fuel_used =  fuelUsed(time2, ff, ff)




for t in time2[:len(time2)-1]:
    cg, mass = CG_time(t, Initial_fuel, ZFM,CG_ZFM,fuel_used)
    cgg.append(cg)
    m.append(mass)
plt.figure()
plt.plot(time2[:len(time2)-1], m)

plt.figure()
plt.plot(time2[:len(time2)-1],cgg)

plt.figure()
plt.plot(time2[:len(time2)-1], list(dict.values(fuel_used)))

plt.show()


