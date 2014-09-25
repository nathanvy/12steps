#/usr/bin/python
#Burgers' Equation in 2D
#
# du/dt + u du/dx + v du/dy = nu (d2u/dx2 + d2u/dy2)
# dv/dt + u dv/dx + v dv/dy = nu (d2v/dx2 + d2v/dy2)
#
# discretized and rearranged:
# u_n+1_i,j = u_n_i,j - u_n_i,j * dt/dx * (u_n_i,j - u_n_i-1,j) - v_n_i,j * dt/dy * (u_n_i,j - u_n_i,j-1) + nu dt/dx^2 * ( u_n_i+1,j - 2 u_n_i,j + u_n_i-1,j) + nu dt/dy^2 * ( u_n_i,j+1 - 2 u_n_i,j + u_n_i,j-1)
# yuck

from mpl_toolkits.mplot3d import Axes3D ## for 3D projected plots
from matplotlib import cm # color map

import numpy
import matplotlib.pyplot

#setup
nx = 41
ny = 41
nt = 120
c = 1
nu = 0.01
dx = 2.0/(nx-1)
dy = 2.0/(ny-1)
sigma = 0.w009
dt = sigma*dx*dy/nu

x = numpy.linspace(0,2,nx)
y = numpy.linspace(0,2,ny)

u = numpy.ones((ny,nx))
un = numpy.ones((ny,nx))
v = numpy.ones((ny,nx))
vn = numpy.ones((ny,nx))
comb = numpy.ones((ny,nx))

u[0.5/dy:1/dy+1, 0.5/dx:1/dx+1] = 2
v[0.5/dy:1/dy+1, 0.5/dx:1/dx+1] = 2

#icfig = matplotlib.pyplot.figure(figsize=(11,7), dpi=100)
#ax = icfig.gca(projection='3d')
#X,Y = numpy.meshgrid(x,y)
#ax.plot_wireframe(X,Y,u[:], cmap=cm.coolwarm)
#ax.plot_wireframe(X,Y,v[:], cmap=cm.coolwarm)
#ax.plot_surface(X,Y,u[:], rstride=1, cstride=1, cmap=cm.coolwarm, linewidth=0, antialiased=False)
#ax.set_xlim(0,2)
#ax.set_ylim(0,2)
#ax.set_zlim(1,2.5)
#ax.zaxis.set_major_locator(LinearLocator(5))
#matplotlib.pyplot.savefig("step8-ic")



#def diffuse(nt):
for n in range(nt+1): ##loop across number of time steps\n",
    un = u.copy()
    vn = v.copy()
    u[1:-1,1:-1] = un[1:-1,1:-1] - dt/dx*un[1:-1,1:-1]*(un[1:-1,1:-1]-un[0:-2,1:-1])-dt/dy*vn[1:-1,1:-1]* (un[1:-1,1:-1]-un[1:-1,0:-2])+nu*dt/dx**2*(un[2:,1:-1]-2*un[1:-1,1:-1]+un[0:-2,1:-1])+ nu*dt/dy**2*(un[1:-1,2:]-2*un[1:-1,1:-1]+un[1:-1,0:-2])
    
    v[1:-1,1:-1] = vn[1:-1,1:-1] - dt/dx*un[1:-1,1:-1]*(vn[1:-1,1:-1]-vn[0:-2,1:-1])-dt/dy*vn[1:-1,1:-1]* \
        (vn[1:-1,1:-1]-vn[1:-1,0:-2])+nu*dt/dx**2*(vn[2:,1:-1]-2*vn[1:-1,1:-1]+vn[0:-2,1:-1])+ \
        nu*dt/dy**2*(vn[1:-1,2:]-2*vn[1:-1,1:-1]+vn[1:-1,0:-2])
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
wire1 = ax.plot_wireframe(X,Y,u[:])
wire2 = ax.plot_wireframe(X,Y,v[:])
#ax.set_zlim(1,2.5)
matplotlib.pyplot.savefig("step8-diffuse-"+'nt')
    
#diffuse(120)
