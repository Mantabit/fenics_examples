import fenics as fe
from fenics import *

#define the constants
T=2.0       #final time
steps=10    #number of time steps
dt=T/steps
alpha=3
beta=1.2

#define the time variable
t=0

#define the mesh and the function space
nx=10
ny=10
mesh=UnitSquareMesh(nx,ny)
V=FunctionSpace(mesh,"P",1)

#define the boundary conditions
def boundary(x,on_boundary):
    return on_boundary
u_D=Expression("1+x[0]*x[0]+alpha*x[1]*x[1]+beta*t",degree=2,alpha=alpha,beta=beta,t=t)
bc=fe.DirichletBC(V,u_D,boundary)
u_n=project(u_D,V)                  #project the inital value of u onto u_n

#define the variational problem
u=TrialFunction(V)
v=TestFunction(V)
f=Constant(beta-2-2*alpha)
a=v*u*dx+dt*dot(grad(u),grad(v))*dx
L=v*(u_n+dt*f)*dx

#compute the solution for each time point. u holds u[n+1] and u_n holds u[n]
u=Function(V)
for n in range(0,steps):
    t+=dt
    u_D.t=t
    
    solve(a==L,u,bc)
    plot(u)
    
    u_n.assign(u)