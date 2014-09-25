#/usr/bin/python
#cavity flow w navier-stokes

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot
import numpy

nx = 41
ny = 41
nt = 500
nit = 50
c = 1
dx = 2.0/(nx-1)
dy = 2.0/(ny-1)
x = numpy.linspace(0,2,nx)
y = numpy.linspace(0,2,ny)
Y,X = numpy.meshgrid(y,x)

rho = 1
nu = 0.1
dt = 0.001

u = numpy.zeros((ny,nx))
v = numpy.zeros((ny,nx))
p = numpy.zeros((ny,nx))
b = numpy.zeros((ny,nx))

#yuck
# this is everything inside the square brackets multiplied by rho in the pressure-poisson equation for cavity flow
def sourceTerms(b, rho, dt, u, v, dx, dy):
    
    b[1:-1,1:-1]=rho*(1/dt*((u[2:,1:-1]-u[0:-2,1:-1])/(2*dx)+(v[1:-1,2:]-v[1:-1,0:-2])/(2*dy))-\
                      ((u[2:,1:-1]-u[0:-2,1:-1])/(2*dx))**2-\
                      2*((u[1:-1,2:]-u[1:-1,0:-2])/(2*dy)*(v[2:,1:-1]-v[0:-2,1:-1])/(2*dx))-\
                      ((v[1:-1,2:]-v[1:-1,0:-2])/(2*dy))**2)
    
    return b
    
def presPoisson(p, dx, dy, b):
    pn = numpy.empty_like(p)
    pn = p.copy()
    
    for q in range(nit):
        pn = p.copy()
        p[1:-1,1:-1] = ((pn[2:,1:-1]+pn[0:-2,1:-1])*dy**2+(pn[1:-1,2:]+pn[1:-1,0:-2])*dx**2)/\
                       (2*(dx**2+dy**2)) -\
                       dx**2*dy**2/(2*(dx**2+dy**2))*b[1:-1,1:-1]
        
        p[-1,:] =p[-2,:]		##dp/dy = 0 at y = 2
        p[0,:] = p[1,:]	 	##dp/dy = 0 at y = 0
        p[:,0]=p[:,1]		   ##dp/dx = 0 at x = 0
        p[:,-1]=0		       ##p = 0 at x = 2
        
    return p

def cavityFlow(nt, u, v, dt, dx, dy, p, rho, nu):
    un = numpy.empty_like(u)
    vn = numpy.empty_like(v)
    b = numpy.zeros((ny, nx))
    
    for n in range(nt):
        un = u.copy()
        vn = v.copy()
        
        b = sourceTerms(b, rho, dt, u, v, dx, dy)
        p = presPoisson(p, dx, dy, b)
        
        u[1:-1,1:-1] = un[1:-1,1:-1]-\
                       un[1:-1,1:-1]*dt/dx*(un[1:-1,1:-1]-un[0:-2,1:-1])-\
                       vn[1:-1,1:-1]*dt/dy*(un[1:-1,1:-1]-un[1:-1,0:-2])-\
                       dt/(2*rho*dx)*(p[2:,1:-1]-p[0:-2,1:-1])+\
                       nu*(dt/dx**2*(un[2:,1:-1]-2*un[1:-1,1:-1]+un[0:-2,1:-1])+\
                           dt/dy**2*(un[1:-1,2:]-2*un[1:-1,1:-1]+un[1:-1,0:-2]))
	
        v[1:-1,1:-1] = vn[1:-1,1:-1]-\
                       un[1:-1,1:-1]*dt/dx*(vn[1:-1,1:-1]-vn[0:-2,1:-1])-\
                       vn[1:-1,1:-1]*dt/dy*(vn[1:-1,1:-1]-vn[1:-1,0:-2])-\
                       dt/(2*rho*dy)*(p[1:-1,2:]-p[1:-1,0:-2])+\
                       nu*(dt/dx**2*(vn[2:,1:-1]-2*vn[1:-1,1:-1]+vn[0:-2,1:-1])+\
                           (dt/dy**2*(vn[1:-1,2:]-2*vn[1:-1,1:-1]+vn[1:-1,0:-2])))
        
        u[0,:] = 0
        u[:,0] = 0
        u[:,-1] = 1
        v[0,:] = 0
        v[-1,:]=0
        v[:,0] = 0
        v[:,-1] = 0
        u[-1,:] = 0
        
    return u, v, p
        
u = numpy.zeros((ny,nx))
v = numpy.zeros((ny,nx))
p = numpy.zeros((ny,nx))
b = numpy.zeros((ny,nx))
nt = 200
u,v,p = cavityFlow(nt, u, v, dt, dx, dy, p, rho, nu)
fig = matplotlib.pyplot.figure(figsize=(11,7),dpi=100)
matplotlib.pyplot.contourf(X,Y,p,alpha=0.5)
matplotlib.pyplot.colorbar()
matplotlib.pyplot.contour(X,Y,p)
matplotlib.pyplot.quiver(X[::2,::2],Y[::2,::2],u[::2,::2],v[::2,::2])
matplotlib.pyplot.xlabel('X')
matplotlib.pyplot.ylabel('Y')
matplotlib.pyplot.savefig('nt200')

u = numpy.zeros((ny,nx))
v = numpy.zeros((ny,nx))
p = numpy.zeros((ny,nx))
b = numpy.zeros((ny,nx))
nt = 700
u,v,p = cavityFlow(nt, u, v, dt, dx, dy, p, rho, nu)
fig = matplotlib.pyplot.figure(figsize=(11,7),dpi=100)
matplotlib.pyplot.contourf(X,Y,p,alpha=0.5)
matplotlib.pyplot.colorbar()
matplotlib.pyplot.contour(X,Y,p)
matplotlib.pyplot.quiver(X[::2,::2],Y[::2,::2],u[::2,::2],v[::2,::2])
matplotlib.pyplot.xlabel('X')
matplotlib.pyplot.ylabel('Y')
matplotlib.pyplot.savefig('nt700')
