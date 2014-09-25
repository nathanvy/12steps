#/usr/bin/python
#2d convection
# du/dt + udu/dx + vdu/dy = 0 and
# dv/dt + udv/dx + vdv/dy = 0
#
# same ICs as for 1D

from mpl_toolkits.mplot3d import Axes3D ## for 3D projected plots
from matplotlib import cm # color map

import numpy
import matplotlib.pyplot

#setup
nx = 101
ny = 101
nt = 80
c = 1
dx = 2.0/(nx-1)
dy = 2.0/(ny-1)
sigma = 0.2
dt = sigma*dx

x = numpy.linspace(0,2,nx)
y = numpy.linspace(0,2,ny)

u = numpy.ones((ny,nx))
v = numpy.ones((ny,nx))
un = numpy.ones((ny,nx))
vn = numpy.ones((ny,nx))

# ICs
u[0.5/dy:1/dy+1, 0.5/dx:1/dx+1] = 2
v[0.5/dy:1/dy+1, 0.5/dx:1/dx+1] = 2

for n in range(nt+1):
    un = u.copy()
    vn = v.copy()

    #u_n+1_i,j = u_n_i,j - u_n_i,j * deltat/deltax * (u_n_i,j - u_n_i-1,j) - v_n_i,j * deltat/deltay * (u_n_i,j - u_n_i,j-1)
    u[1:,1:] = un[1:,1:] - (un[1:,1:]*dt/dx*(un[1:,1:] - un[0:-1,1:])) - vn[1:,1:]*dt/dy*(un[1:,1:] - un[1:,0:-1])
    v[1:,1:] = vn[1:,1:] - (un[1:,1:]*dt/dx*(vn[1:,1:] - vn[0:-1,1:])) - vn[1:,1:]*dt/dy*(vn[1:,1:] - vn[1:,0:-1])

    u[0,:] = 1
    u[-1,:] = 1
    u[:,0] = 1
    u[:,-1] = 1

    v[0,:] = 1
    v[-1,:] = 1
    v[:,0] = 1
    v[:,-1] = 1

fig = matplotlib.pyplot.figure(figsize=(11,7), dpi=100)
ax = fig.gca(projection='3d')
X,Y = numpy.meshgrid(x,y)
ax.plot_surface(X,Y,u, cmap=cm.coolwarm)
matplotlib.pyplot.savefig("step6")
