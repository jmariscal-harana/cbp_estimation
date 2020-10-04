# src-ch1/ODEschemes.py

import numpy as np
from matplotlib.pyplot import plot, show, legend, hold,rcParams,rc, figure, axhline, close,\
    xticks, title, xlabel, ylabel, savefig, axis, grid, subplots, setp

# change some default values to make plots more readable 
LNWDT=3; FNT=10
rcParams['lines.linewidth'] = LNWDT; rcParams['font.size'] = FNT



# define Euler solver
def euler(func, z0, time,q=0,prms=0):
    """The Euler scheme for solution of systems of ODEs.
    z0 is a vector for the initial conditions,
    the right hand side of the system is represented by func which returns
    a vector with the same size as z0. The optionsl source term q is an array (function of time), \
    whereas the optional prms are independent of time  """

    z = np.zeros((np.size(time), np.size(z0)))
    z[0,:] = z0

    for i in range(len(time)-1):
        dt = time[i+1] - time[i]
        z[i+1,:]=z[i,:] + np.asarray(func(z[i,:],time[i],q[i],prms))*dt

    return z

