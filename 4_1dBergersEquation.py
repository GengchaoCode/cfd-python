'''
Solve the Bergers's equation with the finite difference method. The Bergers' equation is a 
nonlinear convection and second-order diffusion equation, just like the Navier-Stokes equation
'''

import numpy as np                              # here we load numpy
import sympy as sym                             # symbolic math in python
from sympy.utilities.lambdify import lambdify as lambdify   # turn symbolic function into a python function
from matplotlib import pyplot as plt            # here we load matplotlib
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator)
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
plt.rcParams["font.family"] = "stix"            # set the font to Times globally
plt.rcParams["mathtext.fontset"] = "stix"       # set the math font to Times


## Define a function to describe the initial condition using sympy
# start by setting up symbolic variables for the three variables in our initial condition
x, nu, t = sym.symbols('x nu t')
phi = sym.exp(-(x-4*t)**2/4/nu/(t+1))+sym.exp(-(x-4*t-2*sym.pi)**2/4/nu/(t+1))
phiprime = phi.diff(x)                          # differentiate phi w.r.t. x

# write the full initial condition and turn it to a python function
uexpres = -2*nu/phi*phiprime+4
ufunc = lambdify((x, nu, t), uexpres)              # the lambda function is a function without a name


## Define the geometry (length), spatial and temporal resolutions
lx = 2                                          # consider a domain that is 2 unit of length
nx = 101                                        # the number of grid points
dx = lx/(nx-1)                                  # grid spacing
nt = 100                                        # the number of time steps
nu = 0.07                                       # fluid kinematic viscosity
dt = dx*nu                                      # time step

# calculate the initial condition based on the input parameters
x = np.linspace(0, 2*np.pi, nx)                 # x-coordinate
u = np.empty(nx)                                # initialize a u for the initial condition
for i in range(nx):
    u[i] = ufunc(x[i], nu, 0)                   # calculate the u at every x


## Plot the initial condition when t = 0
plt.figure(figsize=(3,2))
plt.plot(x, u, marker='.')

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
plt.xlim(0, 2*np.pi)
plt.ylim(-1, 8)
plt.tight_layout(pad=0.1)                       # make the layout tight to minimize the white space

# annotate the current time
ax.annotate('$t = 0.000$ s', xy=(0.75,0.9), xycoords='axes fraction', fontsize=10)

# save and show the figure
folderName = '/home/ygc/Documents/Codes/cfd-python/1dBergers/'
fileName = 'u000.png'
plt.savefig(folderName+fileName, dpi=300)
plt.show(block=False)
plt.pause(0.1)                                  # show the image for 0.1 s
plt.close()


## Solve the non-linear convection and second-order diffusion Bergers' equation using finite difference and plot the results
# finite difference method
un = np.empty(nx)                               # initialize an array to store the current velocity

for t in range(nt):
    un = u.copy()                               # calculate u based on un

    # calculate the value at the interior grid points
    for i in range(1,nx-1):
        u[i] = un[i]-un[i]*dt/dx*(un[i]-un[i-1])+nu*dt/dx**2*(un[i+1]-2*un[i]+un[i-1])

    # set the boundary values according to the periodic boundary condition -> head = end
    u[0] = un[0]-un[0]*dt/dx*(un[0]-un[-2])+nu*dt/dx**2*(un[1]-2*un[0]+un[-2])
    u[-1] = u[0]

    # calculate the analytical solutions
    uright = np.empty(nx)                       # initialize a u for the initial condition
    for i in range(nx):
        uright[i] = ufunc(x[i], nu, (t+1)*dt)   # calculate the u at every x
    
    # plot and compare the numerical and analytical solutions
    plt.figure(figsize=(3,2))
    plt.plot(x, u, marker='.', label='numerical')
    # plt.plot(x, uright, label='analytical')

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
    plt.xlim(0, 2*np.pi)
    plt.ylim(-1, 8)
    plt.tight_layout(pad=0.1)                       # make the layout tight to minimize the white space

    # annotate the current time
    ax.annotate('$t = {0:.3f}$ s'.format((t+1)*dt), xy=(0.75,0.9), xycoords='axes fraction', fontsize=10)
    
    # show the legend
    # plt.legend(fontsize=8)

    # save and show the figure
    fileName = 'u{:0>3d}.png'.format(t+1)
    plt.savefig(folderName+fileName, dpi=300)
    plt.show(block=False)
    plt.pause(0.1)                                  # show the image for 0.1 s
    plt.close()