#from CGcalc import *
import pandas as pd
import numpy as np
from main import data
import matplotlib.pyplot as plt
import matplotlib


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
            momentcg = fuel_moment.iat[i,1] - ((fuel_moment.iat[i,1] - fuel_moment.iat[i-1,1])/(fuel_moment.iat[i,0]-fuel_moment.iat[i-1,0]) )* (fuel_moment.iat[i,0] - fuel)
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
    fuel_mom = interpolatefuel(fuel_mass, fuel_moment)*100
    tot_mass = ZFM+ fuel_mass
    CG = CGshift1(ZFM, CG_ZFM, tot_mass, fuel_mom)
    
    return CG, tot_mass
#-------------------------------
MAC = 80.98*0.0254 # 2.056892[m]
LEMAC = 261.56*0.0254 #  205.6892[m]


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

mom_fuel = interpolatefuel(Initial_fuel, fuel_moment)*100

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
        
    cgg = np.array(cgg)
    
#    f_mom = []
#    for i in range(4800, 5009):
#        f_mom.append(interpolatefuel(i, fuel_moment))
#    
#    
#    plt.figure()
#    plt.plot(range(4800, 5009), f_mom)
#    plt.show()
    
    
#    plt.figure()
#    plt.plot(time[:len(time)-1], m)
#    plt.ylabel('Total aircraft mass [kg]')
#    plt.xlabel('Time [s]')
#    plt.grid()
#    
#    plt.figure()
#    plt.plot(time[:len(time)-1],CG_MAC(cgg*0.0254, LEMAC, MAC))
#    plt.ylabel('Centre of mass [% MAC]')
#    plt.xlabel('Time [s]')
#    plt.grid()
#    
#    plt.show()
    
    ######################################################################
    
    ## Graphing Parameters
    texpsize= [26,28,30]

    SMALL_SIZE  = texpsize[0]
    MEDIUM_SIZE = texpsize[1]
    BIGGER_SIZE = texpsize[2]
    
    plt.style.use('grayscale')
    plt.rc('font', size=MEDIUM_SIZE, family='serif')    ## controls default text sizes
    plt.rc('axes', titlesize=SMALL_SIZE)                ## fontsize of the axes title
    plt.rc('axes', labelsize=SMALL_SIZE)                ## fontsize of the x and y labels
    plt.rc('xtick', labelsize=SMALL_SIZE)               ## fontsize of the tick labels
    plt.rc('ytick', labelsize=SMALL_SIZE)               ## fontsize of the tick labels
    plt.rc('legend', fontsize=SMALL_SIZE)               ## legend fontsize
    plt.rc('figure', titlesize=BIGGER_SIZE)             ## fontsize of the figure title
    plt.rc('text', usetex=False)
    matplotlib.rcParams['lines.linewidth']  = 1.5
    matplotlib.rcParams['figure.facecolor'] = 'white'
    matplotlib.rcParams['axes.facecolor']   = 'white'
    matplotlib.rcParams["legend.fancybox"]  = False
    
    ## Graph
    
    
#    fig, ax = plt.subplots(1,1,squeeze=False,figsize=(16,9))
#    ax[0,0].plot(time[:len(time)-1], m)
#    ax[0,0].set_ylabel(r'Total aircraft mass [kg]')         
#    ax[0,0].set_xlabel(r'Time [s]')
#    ax[0,0].set_xlim(0,time[:len(time)-1][len(time)-2])
#    ax[0,0].grid(True,which="major",color="#999999")
#    ax[0,0].grid(True,which="minor",color="#DDDDDD",ls="--")
#    ax[0,0].minorticks_on()
#    ax[0,0].tick_params(which='major', length=10, width=2, direction='inout')
#    ax[0,0].tick_params(which='minor', length=5, width=2, direction='in')
##    ax[0,0].legend(loc=0, framealpha=1.0).get_frame().set_edgecolor('k')
#    fig.savefig("mass-time-plot2.png", bbox_inches='tight')                                    ## Insert save destination
#    
#    ## If you want to see the figure, else disable last two lines.
#    fig.tight_layout()
#    plt.show()            
#
#
#    fig, ax = plt.subplots(1,1,squeeze=False,figsize=(16,9))
#    ax[0,0].plot(time[:len(time)-1],CG_MAC(cgg*0.0254, LEMAC, MAC))
#    ax[0,0].set_ylabel(r'Centre of mass [% MAC]')   
#    ax[0,0].set_xlabel(r'Time [s]')
#    ax[0,0].set_xlim(0,time[:len(time)-1][len(time)-2])
#    ax[0,0].grid(True,which="major",color="#999999")
#    ax[0,0].grid(True,which="minor",color="#DDDDDD",ls="--")
#    ax[0,0].minorticks_on()
#    ax[0,0].tick_params(which='major', length=10, width=2, direction='inout')
#    ax[0,0].tick_params(which='minor', length=5, width=2, direction='in')
##    ax[0,0].legend(loc=0, framealpha=1.0).get_frame().set_edgecolor('k')
#    fig.savefig("CoM-time-plot2.png", bbox_inches='tight')                                    ## Insert save destination
#    
#    ## If you want to see the figure, else disable last two lines.
#    fig.tight_layout()
#    plt.show()       
    
    cg_full_fuel =  [CG_MAC(7.146208878242001, LEMAC, MAC)]* (len(time)-1)
    cg_zfm = [CG_MAC(CG_ZFM*0.0254, LEMAC, MAC)] * (len(time)-1)
        
    fig, ax = plt.subplots(1,1,squeeze=False,figsize=(16,9))
    ax[0,0].plot(time[:len(time)-1],CG_MAC(cgg*0.0254, LEMAC, MAC))
    ax[0,0].plot(time[:len(time)-1],cg_full_fuel, linewidth = 4, linestyle = 'dashed', color = 'black')
    ax[0,0].plot(time[:len(time)-1],cg_zfm,  linewidth = 4, linestyle = 'dashed', color = 'black')
    
    
    ax[0,0].set_ylabel(r'Centre of mass [% MAC]')   
    ax[0,0].set_xlabel(r'Time [s]')
    ax[0,0].set_xlim(0,time[:len(time)-1][len(time)-2])
    ax[0,0].grid(True,which="major",color="#999999")
    ax[0,0].grid(True,which="minor",color="#DDDDDD",ls="--")
    ax[0,0].minorticks_on()
    ax[0,0].tick_params(which='major', length=10, width=2, direction='inout')
    ax[0,0].tick_params(which='minor', length=5, width=2, direction='in')
    #    ax[0,0].legend(loc=0, framealpha=1.0).get_frame().set_edgecolor('k')
    fig.savefig("CG_range_plot.png", bbox_inches='tight')                                    ## Insert save destination
    
    ## If you want to see the figure, else disable last two lines.
    fig.tight_layout()
    plt.show()  







