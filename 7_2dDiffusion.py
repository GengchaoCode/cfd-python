'''
Solve the 2-D diffusion equation using the finite difference method.
'''

import numpy as np                              # here we load numpy
from matplotlib import pyplot as plt            # here we load matplotlib
from matplotlib import cm                       # colormap
import time, sys                                # here we load some utilities
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator)
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.mplot3d import Axes3D         # new library required for projected 3D plots
plt.rcParams["font.family"] = "stix"            # set the font to Times globally
plt.rcParams["mathtext.fontset"] = "stix"       # set the math font to Times

## Variable declarations
nx = 81                                         # grid points in x-direction
ny = 81                                         # grid points in y-direction
nt = 100                                        # number of time steps
dx = 2/(nx-1)                                   # spatial resolution in x-direction
dy = 2/(ny-1)                                   # spatial resolution in y-direction
nu = 0.05                                       # viscosity
sigma = 0.2                                     # for CFL condition
dt = sigma*dx*dy/nu

x = np.linspace(0,2,nx)                         # x-coordinates
y = np.linspace(0,2,ny)                         # y-coordinates

## Assign initial conditions
# u = 2 when x and y are between 0.5 and 1 and u = 1 everywhere else
u = np.ones((ny,nx))                            # col (x) will always be the last dimension
u[int(0.5/dy):int(1/dy+1), int(0.5/dx):int(1/dx+1)] = 2

## Plot the initial condition
fig = plt.figure()
ax = fig.gca(projection='3d')
X, Y = np.meshgrid(x, y)
surf = ax.plot_surface(X, Y, u, cmap=cm.viridis, antialiased=False)

# set the axis properties
ax.xaxis.set_major_formatter(FormatStrFormatter('%g'))
ax.yaxis.set_major_formatter(FormatStrFormatter('%g'))
ax.zaxis.set_major_formatter(FormatStrFormatter('%g'))
ax.tick_params(labelsize=8)

# set the figure properties
plt.xlabel('$x$ (m)', fontsize=10)
plt.ylabel('$y$ (m)', fontsize=10)
ax.set_zlabel('$u$ (m/s)', fontsize=10)
plt.xlim(0, 2)
plt.ylim(0, 2)
ax.set_zlim(1, 2)
plt.tight_layout(pad=0.1)                       # make the layout tight to minimize the white space

# annotate the current time
ax.annotate('$t = 0.000$ s', xy=(0.75,0.9), xycoords='axes fraction', fontsize=10)

# save and show the figure
folderName = '/home/ygc/Documents/Codes/cfd-python/2dDiffusion/'
fileName = 'u000.png'
plt.savefig(folderName+fileName, dpi=300)
plt.show(block=False)
plt.pause(0.1)                                  # show the image for 0.1 s
plt.close()

## Solve using finite difference and plot the results
for n in range(nt):
    un = u.copy()

    # for i in range(1,nx-1):
    #     for j in range(1,ny-1):
    #         u[j,i] = un[j,i]+nu*dt/dx/dx*(un[j,i+1]-2*un[j,i]+un[j,i-1])+nu*dt/dy/dy*(un[j+1,i]-2*un[j,i]+un[j-1,i])

    # the array operation instead of loops should be much more efficient
    u[1:-1,1:-1] = un[1:-1,1:-1]+nu*dt/dx/dx*(un[1:-1,2:]-2*un[1:-1,1:-1]+un[1:-1,0:-2])+nu*dt/dy/dy*(un[2:,1:-1]-2*un[1:-1,1:-1]+un[0:-2,1:-1])

    # set the boundary values according to the BCs
    u[0,:] = 1
    u[-1,:] = 1
    u[:,0] = 1
    u[:,-1] = 1

    # plot the results
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    surf = ax.plot_surface(X, Y, u, cmap=cm.viridis, antialiased=False)

    # set the axis properties
    ax.xaxis.set_major_formatter(FormatStrFormatter('%g'))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%g'))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%g'))
    ax.tick_params(labelsize=8)

    # set the figure properties
    plt.xlabel('$x$ (m)', fontsize=10)
    plt.ylabel('$y$ (m)', fontsize=10)
    ax.set_zlabel('$u$ (m/s)', fontsize=10)
    plt.xlim(0, 2)
    plt.ylim(0, 2)
    ax.set_zlim(1, 2)
    plt.tight_layout(pad=0.1)                       # make the layout tight to minimize the white space

    # annotate the current time
    ax.annotate('$t = {0:.3f}$ s'.format((n+1)*dt), xy=(0.75,0.9), xycoords='axes fraction', fontsize=10)

    # save and show the figure
    fileName = 'u{:0>3d}.png'.format(n+1)
    plt.savefig(folderName+fileName, dpi=300)
    plt.show(block=False)
    plt.pause(0.1)                                  # show the image for 0.1 s
    plt.close()