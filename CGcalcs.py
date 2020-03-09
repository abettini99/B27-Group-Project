from CGcalc import *
from main import data

time = data['time']
FFr = data['rh_engine_FMF']
FFl = data['lh_engine_FMF']





def CGshift(Mi, CGi, Mf, moment_change):
    CG = (Mi*CGi + moment_change)/Mf
    
    return CG

def CG_MAC(CG_datum):
    CG = (CG_datum - LEMAC)/MAC *100
    return CG


MAC = 2.0569 #[m]
LEMAC = 6.64083 #[m]

# calculate BEM
#--taken from weight measurements
BEM = 9165 #[lbs]
CG_BEM = 291.647954 #inches



#Calculate Payload weight:
#--using passenger/payload data
x_seats = [131,131,214,214,251,251,288,288,170]
w_seats = [180,180,170,187,190,180,178,180,190]

x_bag = [74, 321, 338]
w_bag = [50, 70, 65]

w_PL = sum(w_seats) + sum(w_bag)

mom_PL  = sum([x_seats[i]+w_seats[i] for i in range(len(x_seats))]) + sum([x_bag[j]+w_bag[j] for j in range(len(x_bag))])

#OEW CG:
OEW = BEM+w_PL
CG_OEW = CGshift(BEM, CG_BEM, OEW, mom_PL) 

#calculate ramp weight CG:
#--Using interpolated fuel-moment chart

Initial_fuel = 1000 #[lbs]

mom_fuel = interpolatefuel(Initial_fuel)

RW = OEW + Initial_fuel

CG_RW = CGshift(OEW, CG_OEW, RW, mom_fuel)


#Calculate CG at specific time:
#--Parse fuel-time data
#--calcukate fuel moment using interpolated fuel-moment data
#--Calc new CG




