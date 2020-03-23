# -*- coding: utf-8 -*-
"""
Created on Mon Mar 23 14:29:19 2020

@author: chang
"""

from main import *

def eigenval_symm():
    
    A = -1 * (2 * muc) * (CZadot - 2 * muc) * (2 * muc * KY2**2)
    
    B = ( CXu * (CZadot - 2 * muc) * 2 * muc * KY2**2 
        - 2 * muc * CZa * 2 * muc * KY2 
        + 2 * muc * (CZadot - 2 * muc) * Cmq 
        - (CZq + 2 * muc) * Cmadot * 2 * muc )
    
    C = ( CXu * CZa * (2 * muc * KY2) 
         - CXu * (CZadot - 2 * muc) * Cmq 
         + 2 * muc * CZa * Cmq 
         - CZu * Cmadot * CXq 
        + CXq * (CZadot - 2 * muc) * Cmu 
        - (CZq + 2 * muc) * Cma * 2 * muc 
        + (CZq + 2 * muc) * Cmadot * CXu 
        - 2 * muc * KY2**2 * CXa * CZu
        + CX0 * Cmadot * 2 * muc )
    
    D = ( -1 * CXu * CZa * Cmq 
         - CZu * Cma * CXq 
         - Cmu * CXa * (CZq + 2 * muc) 
         + CXq * CZa * Cmu
        + (CZq + 2 * muc) * Cma * CXu 
        + Cmq * CXa * CZu 
        - CZu * Cmadot * CZ0 
        + CZ0 * (CZadot - 2 * muc) * Cmu
        + CX0 * Cma * 2 * muc 
        - CX0 * Cmadot * CXu )
    
    E = ( -1 * CZu * Cma * CZ0 
        + Cmu * CXa * CX0 
        + CZ0 * CZa * Cmu 
        - CX0 * Cma * CXu )
    
    return np.roots([A,B,C,D,E])
    

def eigenval_phugoid():
    
    A = ( 2 * muc * CZa * Cmq 
        - (CZq + 2 * muc) * Cma * 2 * muc )
    
    B = ( -1 * CXu * CZa * Cmq 
         - CZu * Cma * CXq 
         - Cmu * CXa * (CZq + 2 * muc) 
         + CXq * CZa * Cmu
         + (CZq + 2 * muc) * Cma * CXu 
         + Cmq * CXa * CZu 
         + CX0 * Cma * 2 * muc )
    
    C = ( -1 * CZu * Cma * CZ0 
         + Cmu * CXa * CX0 
         + CZ0 * CZa * Cmu 
         - CX0 * Cma * CXu )
    
    return np.roots([A,B,C])


def eigenval_shortperiod():
    
    A = -1 * (CZadot - 2 * muc) * 2 * muc * KY2**2
    
    B = ( -1 * CZa * 2 * muc * KY2**2
         + (CZadot - 2 * muc) * Cmq
         + (CZq + 2 * muc) * Cmadot )
    
    C = ( -1 * CZa * Cmq 
         - CX0 * Cmadot
         + (CZq + 2 * muc) * Cma )
    
    D = -1 * CX0 * Cma
    
    return np.roots([A,B,C,D])


def eigenval_assym():
    
    A = 0.5 * ( -1 * (CYbdot - 2 * mub) * 4 * mub * KX2**2 * 4 * mub * KZ2**2
               + (4 * mub * KXZ)**2 * (CYbdot - 2 * mub) )
    
    B = 0.5 * ( -1 * CYb * 4 * mub * KX2**2 * 4 * mub * KZZ**2 
               + (CYbdot - 2 * mub) * Clp * 4 * mub  * KZ2**2
               + (CYbdot - 2 * mub) * 4 * mub * KX2 **2 * Cnr
               - Cnbdot * CYp * 4 * mub * KXZ
               - (CYr - 4 * mub) * 4 * mub * KX2**2 * Cnbdot
               + Clr * 4 * mub * KXZ * (CYbdot - 2 * mub)
               + 4 * mub * KXZ * Cnp * (CYbdot - 2 * mub)
               + (4 * mub * KXZ)**2 * CYb )
    
    C = ( 0.5 * ( CYb * Clp * 4 * mub * KZ2 ** 2
            + CYb * 4 * mub * KX2**2 * Cnr
            - (CYbdot - 2 * mub) * Clp * Cnr
            - Clb * 4 * mub * KXZ * (CYr - 4 * mub)
            - Cnb * CYp * 4 * mub * KXZ
            - Cnbdot * CYp * Clr
            + (CYr - 4 * mub) * Clp * Cnbdot
            - (CYr - 4 * mub) * 4 * mub * KX2**2 * Cnb
            + Clr * Cnp * (CYbdot - 2 * mub)
            + Clr * 4 * mub * KXZ * CYb
            + 4 * mub * KXZ * Cnp * CYb
            - 4 * mub * KZ2**2 * CYp * Clb )
         - Cnbdot * CL * 4 * mub * KXZ )
    
    D = ( -0.5 * ( -1 * CYb * Clp * Cnr
                - Clb * Cnp * (CYr - 4 * mub)
                - Cnb * CYp * Clr
                + (CYr - 4 * mub) * Clp * Cnb
                + Clr * Cnp * CYb
                + Cnr * CYp * Clb)
        - ( Cnb * CL * 4 * mub * KXZ 
           + Cnbdot * CL * Clr
           + 4 * mub * KZ2**2 * CL * Clb ) )
    
    E = (-1 * Cnb * CL * Clr
         + Cnr * CL * Clb)
    
    return np.roots([A,B,C,D,E])


