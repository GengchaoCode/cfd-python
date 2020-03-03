'''
Solve the 2-D Cavity Flow problem using the finite difference method.
'''
import os
os.system('clear')                                  # clear the screen
import numpy as np                                  # here we load numpy
from matplotlib import pyplot as plt                # here we load matplotlib
from matplotlib import cm                           # colormap
import time, sys                                    # here we load some utilities
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator)
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.mplot3d import Axes3D             # new library required for projected 3D plots
plt.rcParams["font.family"] = "Times New Roman"     # set the font to Times globally
plt.rcParams["mathtext.fontset"] = "stix"           # set the math font to Times

#==================================#
#----- Self-defined Functions -----#
#==================================#
# Construct a function to solve the pressure Poisson equation
def pressure_poisson(p, u, v, rho, dt, dx, dy, nnt):
    # iterate until the pressure field fulfill the governing equation
    for nn in range(nnt):
        pn = p.copy()

        a = dy**2*(pn[1:-1,2:]+pn[1:-1,:-2])+dx**2*(pn[2:,1:-1]+pn[:-2,1:-1])/2/(dx**2+dy**2)
        b = rho*(dx**2)*(dy**2)/2/(dx**2+dy**2)
        c = 1/dt*((u[1:-1,2:]-u[1:-1,:-2])/(2*dx)+(v[2:,1:-1]-v[:-2,1:-1])/(2*dy))
        d = ((u[1:-1,2:]-u[1:-1,:-2])/2/dx)**2
        e = 2*((u[2:,1:-1]-u[:-2,1:-1])/2/dy)*((v[1:-1,2:]-v[1:-1,:-2])/2/dx)
        f = ((v[2:,1:-1]-v[:-2,1:-1])/2/dy)**2

        p[1:-1,1:-1] = a+b*(-c+d+e+f)

    # set the boundary values according to the boundary condition
    p[:,0] = p[:,1]                                 # dp/dx = 0 @ x = 0
    p[:,-1] = p[:,-2]                               # dp/dx = 0 @ x = 2
    p[-1,:] = 0                                     # p = 0 @ y = 2 (free surface)
    p[0,:] = p[1,:]                                 # dp/dy = 0 @ y = 0

    return p                                        # return the updated pressure field

# Construct a function to solve the momentum equation
def ns_solver(p, u, v, rho, dt, dx, dy, nu):
    un = u.copy()
    vn = v.copy()

    # solve the horizontal velocity at the interior points using FDM
    ua = dt/dx*un[1:-1,1:-1]*(un[1:-1,1:-1]-un[1:-1,:-2])
    ub = dt/dy*vn[1:-1,1:-1]*(un[1:-1,1:-1]-un[:-2,1:-1])
    uc = dt/2/rho/dx*(p[1:-1,2:]-p[1:-1,:-2])
    ud = nu*dt/dx/dx*(un[1:-1,2:]-2*un[1:-1,1:-1]+un[1:-1,:-2])
    ue = nu*dt/dy/dy*(un[2:,1:-1]-2*un[1:-1,1:-1]+un[:-2,1:-1])

    u[1:-1,1:-1] = un[1:-1,1:-1]-ua-ub-uc+ud+ue

    # set the boundary values for the horizontal velocity according to the boundary condition
    u[:,0] = 0                                      # u = 0 @ x = min
    u[:,-1] = 0                                     # u = 0 @ x = max
    u[0,:] = 0                                      # u = 0 @ y = min
    u[-1,:] = 1                                     # u = 1 @ y = max (the lid)

    # solve the vertical velocity at the interior points using FDM
    va = dt/dx*un[1:-1,1:-1]*(vn[1:-1,1:-1]-vn[1:-1,:-2])
    vb = dt/dy*vn[1:-1,1:-1]*(vn[1:-1,1:-1]-vn[:-2,1:-1])
    vc = dt/2/rho/dy*(p[2:,1:-1]-p[:-2,1:-1])
    vd = nu*dt/dx/dx*(vn[1:-1,2:]-2*vn[1:-1,1:-1]+vn[1:-1,:-2])
    ve = nu*dt/dy/dy*(vn[2:,1:-1]-2*vn[1:-1,1:-1]+vn[:-2,1:-1])

    v[1:-1,1:-1] = vn[1:-1,1:-1]-va-vb-vc+vd+ve

    # set the boundary values for the vertical velocity according to the boundary condition
    v[:,0] = 0                                      # v = 0 @ x = min
    v[:,-1] = 0                                     # v = 0 @ x = max
    v[0,:] = 0                                      # v = 0 @ y = min
    v[-1,:] = 0                                     # v = 0 @ y = max
    
    return u, v

#=========================#
#----- Major Routine -----#
#=========================#
# Parameter declarations
nx = 41
ny = 41
nt = 600                                            # time steps to envolve the N-S equation
nnt = 50                                            # substeps to iterate the pressure Poisson equation
dx = 2/(nx-1)
dy = 2/(ny-1)
x = np.linspace(0, 2, nx)
y = np.linspace(0, 2, ny)
X, Y = np.meshgrid(x, y)

rho = 1.0
nu = 0.1
dt = 0.002

u = np.zeros((ny,nx))
v = np.zeros((ny,nx))
p = np.zeros((ny,nx))

# Solve the ns equation
for n in range(nt):
    p = pressure_poisson(p, u, v, rho, dt, dx, dy, nnt)
    u, v = ns_solver(p, u, v, rho, dt, dx, dy, nu)

    #===========================#
    #----- Post-processing -----#
    #===========================#
    # Plot the results
    plt.figure(figsize=(4,3.4))
    plt.contourf(X, Y, p, levels=np.linspace(-2.4,2.4,9), alpha=0.5, cmap=cm.viridis)   # plot the pressure field as a contour

    # Set the colorbar properties
    cbar = plt.colorbar(fraction=0.05, pad=0.02)        # show the colorbar
    cax = cbar.ax                                       # get the axis of the color bar
    cax.set_ylabel('Pressure, $p$', fontsize=10)
    # cax.yaxis.set_major_formatter(FormatStrFormatter('%g'))
    cax.tick_params(axis='y', which='major', direction='in', labelsize=8)

    plt.contour(X, Y, p, cmap=cm.viridis)               # show the contour lines
    plt.quiver(X[::2,::2], Y[::2,::2], u[::2,::2], v[::2,::2])  # plot the velocity field

    # Set the axis properties
    ax = plt.gca()
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.xaxis.set_major_formatter(FormatStrFormatter('%g'))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%g'))
    ax.tick_params(axis='both', which='major', direction='in', labelsize=8)
    ax.tick_params(axis='both', which='minor', direction='in')

    # Set the figure property
    plt.title('Driven Cavity Flow'.format(fontsize=12))
    plt.xlabel('$x$', fontsize=10)
    plt.ylabel('$y$', fontsize=10)
    plt.axis('scaled')                                  # make the axis equal
    plt.xlim(0, 2)
    plt.ylim(0, 2)
    plt.xticks(np.arange(0,2.1,0.5))
    plt.yticks(np.arange(0,2.1,0.5))
    plt.tight_layout(pad=0.1)                           # make the layout tight to minimize the white space

    # Save and close the figure
    folderName = './2dCavity/'
    if not os.path.isdir(folderName):
        os.mkdir(folderName)

    figName = 'cavity_{:0>3d}.png'.format(n+1)
    print('Saving figure '+figName)
    plt.savefig(folderName+figName, dpi=300)            # save the figure
    plt.close()