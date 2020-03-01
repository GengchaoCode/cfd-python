'''
Solve the 2-D Laplace equation using the finite difference method. In stead of tracking a wave through time,
the Laplace equation calculates the equilibrium state of a system under the supplied boundary condition.
'''

import numpy as np                                  # here we load numpy
from matplotlib import pyplot as plt                # here we load matplotlib
from matplotlib import cm                           # colormap
import time, sys                                    # here we load some utilities
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator)
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.mplot3d import Axes3D             # new library required for projected 3D plots
plt.rcParams["font.family"] = "Time New Roman"      # set the font to Times globally
plt.rcParams["mathtext.fontset"] = "stix"           # set the math font to Times

## Define a function to plot the 2D results in a 3D figure
def plot2d(x, y, p, t):
    """A function to present the 2-D simulation results in a 3-D way.
    
    Arguments:
        x {meshgrid} -- x coordinates
        y {meshgrid} -- y coordinates
        p {meshgrid} -- pressure value
        t {float} -- time
    """
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    surf = ax.plot_surface(x, y, p, cmap=cm.viridis, antialiased=False)

    # set the axis properties
    ax.xaxis.set_major_formatter(FormatStrFormatter('%g'))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%g'))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%g'))
    ax.tick_params(labelsize=8)

    # set the figure properties
    plt.xlabel('$x$ (m)', fontsize=10)
    plt.ylabel('$y$ (m)', fontsize=10)
    ax.set_zlabel('$p$ (Pa)', fontsize=10)
    plt.xlim(0, 2)
    plt.ylim(0, 1)
    # ax.set_zlim(1, 2)
    plt.tight_layout(pad=0.1)                       # make the layout tight to minimize the white space

    # annotate the current time
    ax.annotate('$t = {0:.3f}$ s'.format(t), xy=(0.75,0.9), xycoords='axes fraction', fontsize=10)

    # show and save the figure
    plt.show()

## Define a function to solve the Laplace equation using finite difference
def laplace2d(p, dx, dy, diffTarget):
    """Iterate the finite difference form of the Laplace equation.
    
    Arguments:
        p {2D array} -- pressure field
        dx {float} -- space step in x
        dy {float} -- space step in y
        diffTarget {float} -- threshold for convergence
    """
    diff = 1
    pn = np.empty_like(p)

    while diff>diffTarget:
        pn = p.copy()
        p[1:-1,1:-1] = (dy**2*(pn[1:-1,2:]+pn[1:-1,:-2]) \
                      + dx**2*(pn[2:,1:-1]+pn[:-2,1:-1])) \
                      / (2*(dx**2+dy**2))
        
        # set the boundary values - done during the calculation
        p[:,0] = 0                                  # p = 0 @ x = 0
        p[:,-1] = np.arange(0,1,dy)                 # p = y @ x = 2
        p[0,:] = p[1,:]                             # dp/dy = 0 @ y = 0
        p[-1,:] = p[-2,:]                           # dp/dy = 0 @ y = 1

        # update the difference two consective time steps based on L1 Norm
        diff = np.sum(np.abs(p)-np.abs(pn))/np.sum(np.abs(pn))

## Variable declarations
nx = 41                                             # nodes in x-direction
ny = 41                                             # nodes in y-direction
lx = 2                                              # domain extent in x
ly = 1                                              # domain extent in y
dx = lx/(nx-1)                                      # space step in x
dy = ly/(ny-1)                                      # space step in y

## Define the initial condition (boundary condition is applied during calculation)
p = np.zeros((ny,nx))                               # create a x by y array of 0s
x = np.linspace(0, lx, nx)
y = np.linspace(0, ly, ny)
p[:,-1] = y

plot2d(x, y, p, 0)
