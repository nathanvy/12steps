#/usr/bin/python
#2D Poisson Eqn
#
#d2p/dx2 + d2p/dy2 = b
#
#discretize
# (p_n_i+1,j - 2p_n_i,j + p_n_i-1,j)/deltaX^2 + (p_n_i,j+1 - 2 p_n_i,j + p_n_i,j-1)/deltaT^2 = b_n_i,j
#
#solve for p_n_i,j
# p_n_i,j = [ (p_n_i+1,j + p_n_i-1,j)(deltaY^2) + (p_n_i,j+1 + p_n_i,j-1)(deltaX^2) - b_n_i,j deltaX^2 deltaY^2 ]/ 2 (deltaX^2 + deltaY^2)
#
#IC p=0 everywhere
#BC p=0 at x=0,2 and y=0,1
#source term b_ij = 100 at i=nx/4, j=ny/4
#b_ij = -100 at i= 3nx/4, j=3ny/4
#b_ij = 0 elsewhere

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot
import numpy

def plot2d(x, y, p, name):
    fig = matplotlib.pyplot.figure(figsize=(11,7), dpi=100)
    ax = fig.gca(projection='3d')
    X,Y = numpy.meshgrid(x,y)
    surf = ax.plot_surface(X, Y, p[:], rstride=1, cstride=1, cmap=cm.coolwarm, linewidth=0, antialiased=False)
    ax.set_xlim(0,2)
    ax.set_ylim(0,1)
    ax.view_init(30,225)
    matplotlib.pyplot.savefig(name)


nx = 50
ny = 50
nt = 100
xmin = 0.
xmax = 2.
ymin = 0.
ymax = 1.

dx = (xmax-xmin)/(nx-1)
dy = (ymax-ymin)/(ny-1)

#init
p = numpy.zeros((nx, ny))
pd = numpy.zeros((nx, ny))
b = numpy.zeros((nx, ny))
x = numpy.linspace(xmin, xmax, nx)
y = numpy.linspace(ymin, ymax, ny)

# source
b[nx/4][ny/4] = 100
b[3*nx/4][3*ny/4] = -100

plot2d(x,y,p,"poisson-ic")

for it in range(nt):
    pd[:][:] = p[:][:]

    p[1:nx-1,1:ny-1] = ( dy**2/(2*(dx**2+dy**2))*(pd[2:nx,1:ny-1]+pd[0:nx-2,1:ny-1]) + 
                         dx**2/(2*(dx**2+dy**2))*(pd[1:nx-1,2:ny]+pd[1:nx-1,0:ny-2]) -
                         b[1:nx-1,1:ny-1]*dx**2*dy**2/(2*(dx**2+dy**2)) )
    p[0,:] = p[nx-1,:] = p[:,0] = p[:,ny-1] = 0.0

plot2d(x,y,p,"poisson-post")
