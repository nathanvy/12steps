#/usr/bin/python
#2D Laplace eqn
#
#IC: p=0 everywhere
#BC: p=0 at x=0, p=y at x=2, dp/dy=0 at y=0

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

import numpy
import matplotlib.pyplot

def plot2d(x, y, p, name):
    fig = matplotlib.pyplot.figure(figsize=(11,7), dpi=100)
    ax = fig.gca(projection='3d')
    X,Y = numpy.meshgrid(x,y)
    surf = ax.plot_surface(X, Y, p[:], rstride=1, cstride=1, cmap=cm.coolwarm, linewidth=0, antialiased=False)
    ax.set_xlim(0,2)
    ax.set_ylim(0,1)
    ax.view_init(30,225)
    matplotlib.pyplot.savefig(name)
    
def laplace2d(p, y, dx, dy, l1norm_target):
    l1norm=1
    pn = numpy.empty_like(p)
    
    while l1norm > l1norm_target:
        pn = p.copy()
        p[1:-1,1:-1] = (dy**2*(pn[2:,1:-1]+pn[0:-2,1:-1])+dx**2*(pn[1:-1,2:]+pn[1:-1,0:-2]))/(2*(dx**2+dy**2))
        p[0,0] = (dy**2*(pn[1,0]+pn[-1,0])+dx**2*(pn[0,1]+pn[0,-1]))/(2*(dx**2+dy**2))
        p[-1,-1] = (dy**2*(pn[0,-1]+pn[-2,-1])+dx**2*(pn[-1,0]+pn[-1,-2]))/(2*(dx**2+dy**2))
        p[:,0] = 0        ##p = 0 at x = 0
        p[:,-1] = y       ##p = y at x = 2
        p[0,:] = p[1,:]   ##dp/dy = 0 at y = 0
        p[-1,:] = p[-2,:] ##dp/dy = 0 at y = 1

        l1norm = (numpy.sum(numpy.abs(p[:])-numpy.abs(pn[:])))/numpy.sum(numpy.abs(pn[:]))
        
    return p

nx = 31
ny = 31
c = 10
dx = 2.0/(nx-1)
dy = 2.0/(ny-1)

#ic's
p = numpy.zeros((ny,nx))

#plotting aids
x = numpy.linspace(0,2,nx)
y = numpy.linspace(0,1,ny)

#bc's
p[:,0] = 0        ##p = 0 at x = 0
p[:,-1] = y       ##p = y at x = 2
p[0,:] = p[1,:]   ##dp/dy = 0 at y = 0
p[-1,:] = p[-2,:] ##dp/dy = 0 at y = 1

plot2d(x,y,p,"ICplot")
laplace2d(p, y, dx, dy, .01)
plot2d(x,y,p,"l1target")
