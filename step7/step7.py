#/usr/bin/python
#2D diffusion
#
# du/dt = nu d2u/dx2 + nu d2u/dy2
#discretizes to:
#
# (u_n+1_i,j - u_n_i,j)/deltaT = nu( u_n_i+1,j - 2u_n_i,j + u_n_i-1,j)/deltaX + nu( u_n_i,j+1 - 2u_n_i,j + u_n_i,j-1)
# solving for the only unknown:
# u_n+1_i,j = u_n_i,j + nu dt/dx (u_n_i+1,j - 2u_n_i,j + u_n_i-1,j) + nu dt/dy(u_n_i,j+1 -2u_n_i,j + u_n_i,j-1)
#

from mpl_toolkits.mplot3d import Axes3D ## for 3D projected plots
from matplotlib import cm # color map

import numpy
import matplotlib.pyplot

#setup
nx = 31
ny = 31
nt = 17
nu = 0.05
dx = 2.0/(nx-1)
dy = 2.0/(ny-1)
sigma = 0.25
dt = sigma*dx

x = numpy.linspace(0,2,nx)
y = numpy.linspace(0,2,ny)

u = numpy.ones((ny,nx))
un = numpy.ones((ny,nx))

u[0.5/dy:1/dy+1, 0.5/dx:1/dx+1] = 2

icfig = matplotlib.pyplot.figure()
ax = icfig.gca(projection='3d')
X,Y = numpy.meshgrid(x,y)
ax.plot_surface(X,Y,u[:], rstride=1, cstride=1, cmap=cm.coolwarm, linewidth=0, antialiased=False)
ax.set_xlim(0,2)
ax.set_ylim(0,2)
ax.set_zlim(1,2.5)
#ax.zaxis.set_major_locator(LinearLocator(5))
matplotlib.pyplot.savefig("step7-ic")

def diffuse(nt):
    u[0.5/dy:1/dy+1, 0.5/dx:1/dx+1] = 2
    for n in range(nt+1):
        un = u.copy()
        u[1:-1,1:-1] = un[1:-1,1:-1] + nu * dt/dx**2 * ( un[2:,1:-1]-2*un[1:-1,1:-1]+un[0:-2,1:-1])+nu*dt/dy**2*(un[1:-1,2:]-2*un[1:-1,1:-1]+un[1:-1,0:-2])
        u[0,:] = 1
        u[:,0] = 1
        u[-1,:] = 1
        u[:,-1] = 1

    fig = matplotlib.pyplot.figure()
    ax = fig.gca(projection='3d')
    surf = ax.plot_surface(X,Y,u[:], rstride=1, cstride=1, cmap=cm.coolwarm, linewidth=0, antialiased=True)
    ax.set_zlim(1,2.5)
    matplotlib.pyplot.savefig("step7-diffuse-"+'nt')

diffuse(10)
