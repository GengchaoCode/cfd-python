'''
Solve the 2-D Laplace equation using the finite difference method. In stead of tracking a wave through time,
the Laplace equation calculates the equilibrium state of a system under the supplied boundary condition.
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

## Define a function to plot the 2D results in a 3D figure
def plot2d(X, Y, p, t):
    """A function to present the 2-D simulation results in a 3-D way.
    
    Arguments:
        X {meshgrid} -- x coordinates
        Y {meshgrid} -- y coordinates
        p {meshgrid} -- pressure values
        t {int} -- iterations
    """
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    surf = ax.plot_surface(X, Y, p, cmap=cm.viridis, antialiased=False)

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
    ax.set_zlim(-0.05, 0.05)
    ax.view_init(30, 225)
    plt.tight_layout(pad=0.1)                       # make the layout tight to minimize the white space

    # annotate the current time
    ax.annotate('Iterations: {}'.format(t+1), xy=(0.75,0.9), xycoords='axes fraction', fontsize=10)

## Parameter Declaration
nx = 51
ny = 51
nt = 200
xmin = 0
xmax = 2
ymin = 0
ymax = 1
dx = (xmax-xmin)/(nx-1)
dy = (ymax-ymin)/(ny-1)

# make the folder for saving figures
folderName = './2dPoisson/'
if not os.path.isdir(folderName):
    os.mkdir(folderName)

## Initialization
p = np.zeros((ny,nx))
pd = np.zeros((ny,nx))

x = np.linspace(xmin, xmax, nx)
y = np.linspace(ymin, ymax, ny)
X, Y = np.meshgrid(x, y)

b = np.zeros((ny,nx))
b[int(ny/4),int(nx/4)] = 100
b[int(3*ny/4),int(3*nx/4)] = -100

## Advance in pseudo-time 
for t in range(nt):
    pd = p.copy()

    p[1:-1,1:-1] = ((pd[1:-1,2:]+pd[1:-1,:-2])*dy**2+(pd[2:,1:-1]+pd[:-2,1:-1])*dx**2-b[1:-1,1:-1]*dx**2*dy**2)/2/(dx**2+dy**2)

    # plot and save the results
    plot2d(X, Y, p, t)
    fileName = 'p{:0>3d}.png'.format(t)
    print('Saving figures: {}'.format(fileName))
    plt.savefig(folderName+fileName, dpi=300)
    plt.close()