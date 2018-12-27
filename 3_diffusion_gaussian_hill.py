import fenics as fe
from fenics import *

#define constants
steps=30
beta=0.5
T=1
dt=T/steps
a=5

#define the time variable
t=0

#define the mesh and trial/test space
nx=ny=30
mesh=RectangleMesh(Point(-2,-2),Point(2,2),nx,ny)
V=FunctionSpace(mesh,"P",1)

#define the boundary and boundary condition
def boundary(x,on_boundary):
    return on_boundary
u_D=Constant(0)
bc=DirichletBC(V,u_D,boundary)

#define the initial condition u(x,y,t=0) and the variable u_n
initial=Expression("exp(-a*x[0]*x[0]-a*x[1]*x[1])",degree=2,a=a)
u_n=interpolate(initial,V)

#create a pvd/vtu file descriptor
vtkfile=File("heat_gaussian/solution.pvd")

#define the variational problem
u=TrialFunction(V)
v=TestFunction(V)
f=Constant(0)
a=u*v*dx+beta*dt*dot(grad(u),grad(v))*dx
L=u_n*v*dx+beta*dt*f*v*dx

u=Function(V)
for n in range(0,steps):
    t+=dt
    
    solve(a==L,u,bc)
    
    vtkfile<<(u,t)
    plot(u)
    plot(mesh)
    
    u_n.assign(u)