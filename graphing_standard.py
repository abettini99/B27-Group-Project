## Library Imports
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

## Define text sizes for **SAVED** pictures (texpsize -- text export size)
texpsize= [26,28,30]

## Input Arrays
x = np.linspace(1,10,10)
y = np.ones((10))*2

## Graphing Parameters
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
fig, ax = plt.subplots(1,1,squeeze=False,figsize=(16,9))
ax[0,0].plot(x, y, label="test")
ax[0,0].plot(x+x, y+y, label="test2", linestyle="dashed")
#ax[0,0].loglog(x, y, marker = "s", color='black', markerfacecolor='none', markeredgewidth=2, markersize=6, label="test")
ax[0,0].set_ylabel(r"Deflection in $y^{\prime}$ direction $u\,\,[m]$")          ## String is treatable as latex code
ax[0,0].set_xlabel(r"Position along Aileron $x\,\,[m]$")
#ax[0,0].set_xlim(0,x[-1])
ax[0,0].grid(True,which="major",color="#999999")
ax[0,0].grid(True,which="minor",color="#DDDDDD",ls="--")
ax[0,0].minorticks_on()
ax[0,0].tick_params(which='major', length=10, width=2, direction='inout')
ax[0,0].tick_params(which='minor', length=5, width=2, direction='in')
ax[0,0].legend(loc=0, framealpha=1.0).get_frame().set_edgecolor('k')
fig.savefig("yeet.png", bbox_inches='tight')                                    ## Insert save destination

## If you want to see the figure, else disable last two lines.
fig.tight_layout()
plt.show()            
