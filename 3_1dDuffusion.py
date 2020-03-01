'''
Solve the 1-D diffusion equation using the finite difference method.
'''

import numpy as np                              # here we load numpy
from matplotlib import pyplot as plt            # here we load matplotlib
import time, sys                                # here we load some utilities
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator)
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
plt.rcParams["font.family"] = "stix"            # set the font to Times globally
plt.rcParams["mathtext.fontset"] = "stix"       # set the math font to Times

## Define the geometry (length), spatial and temporal resolution, and the wavespeed
lx = 2                                          # consider a domain that is 2 unit of length
nx = 41                                         # the number of grid points
dx = lx/(nx-1)                                  # grid spacing
nt = 100                                        # the number of time steps
nu = 0.3                                        # fluid kinematic viscosity
sigma = 0.2                                     # a parameter we will learn later
dt = sigma*dx**2/nu                             # time step for stability

## Define the initial condition - a hat function
u = np.ones(nx)                                 # create an array which is nx elements long with every value equal to 1
ic1 = int(0.5/dx)                               # left bound index
ic2 = int(1/dx)+1                               # right bound index
u[ic1:ic2] = 2                                  # u0 = 2 when x is between 0.5 and 1 and u0 = 1 everywhear else

## Plot the initial condition when t = 0
plt.figure(figsize=(3,2))
plt.plot(np.linspace(0,2,nx), u)

# set the axis properties
ax = plt.gca()
ax.xaxis.set_minor_locator(AutoMinorLocator(2))
ax.yaxis.set_minor_locator(AutoMinorLocator(2))
ax.xaxis.set_major_formatter(FormatStrFormatter('%g'))
ax.yaxis.set_major_formatter(FormatStrFormatter('%g'))
ax.tick_params(axis='both', which='major', direction='in', labelsize=8)
ax.tick_params(axis='both', which='minor', direction='in')

# set the figure properties
plt.xlabel('$x$ (m)', fontsize=10)
plt.ylabel('$u$ (m/s)', fontsize=10)
plt.xlim(0, 2)
plt.ylim(0.9, 2.1)
plt.tight_layout(pad=0.1)                       # make the layout tight to minimize the white space

# annotate the current time
ax.annotate('$t = 0.000$ s', xy=(0.75,0.9), xycoords='axes fraction', fontsize=10)

# save and show the figure
folderName = '/home/ygc/Documents/Codes/cfd-python/1dDiffusion/'
fileName = 'u000.png'
plt.savefig(folderName+fileName, dpi=300)
plt.show(block=False)
plt.pause(0.1)                                  # show the image for 0.1 s
plt.close()

## Solve the convection equation using the finite difference method and plot the result
un = np.ones(nx)                                # initialize an array to store the current velocity

for n in range(nt):                             # advance nt cycles in time
    un = u.copy()                               # copy the current velocity into un
    for i in range(1,nx-1):                       # the left-most value is not updated
        u[i] = un[i]+nu*dt/dx/dx*(un[i+1]-2*un[i]+un[i-1])

    # Plot the velocity profile
    plt.figure(figsize=(3,2))
    plt.plot(np.linspace(0,2,nx), u)

    # set the axis properties
    ax = plt.gca()
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.xaxis.set_major_formatter(FormatStrFormatter('%g'))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%g'))
    ax.tick_params(axis='both', which='major', direction='in', labelsize=8)
    ax.tick_params(axis='both', which='minor', direction='in')

    # set the figure properties
    plt.xlabel('$x$ (m)', fontsize=10)
    plt.ylabel('$u$ (m/s)', fontsize=10)
    plt.xlim(0, 2)
    plt.ylim(0.9, 2.1)
    plt.tight_layout(pad=0.1)                       # make the layout tight to minimize the white space

    # annotate the current time
    ax.annotate('$t = {0:.3f}$ s'.format((n+1)*dt), xy=(0.75,0.9), xycoords='axes fraction', fontsize=10)

    # save and show the figure
    fileName = 'u{:0>3d}.png'.format(n+1)
    plt.savefig(folderName+fileName, dpi=300)
    plt.show(block=False)
    plt.pause(0.1)                                  # show the image for 0.1 s
    plt.close()