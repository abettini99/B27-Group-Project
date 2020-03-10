from CGcalcs import *

MAC = 2.0569 #[m]
LEMAC = 6.64083 #[m]


#--Basic empty mass (BEM) taken from weight measurements
BEM = 9165 #[lbs]
CG_BEM = 291.647954 #inches


#--Payload weight using passenger/payload data
x_seats = [131,131,214,214,251,251,288,288,170]
M_seats = [180,180,170,187,190,180,178,180,190]

x_bag = [74, 321, 338]
M_bag = [50, 70, 65]

M_PL = sum(M_seats) + sum(M_bag)

mom_PL  = sum([x_seats[i]*M_seats[i] for i in range(len(x_seats))]) + sum([x_bag[j]*M_bag[j] for j in range(len(x_bag))])

#ZFW CG:
ZFM = BEM+M_PL
CG_ZFM = CGshift(BEM, CG_BEM, ZFM, mom_PL) 

#calculate ramp weight CG:
#--Using interpolated fuel-moment chart

Initial_fuel = 2000 #[lbs]

mom_fuel = interpolatefuel(Initial_fuel, fuel_moment)

RM = ZFM + Initial_fuel

CG_RM = CGshift(ZFM, CG_ZFM, RM, mom_fuel)
