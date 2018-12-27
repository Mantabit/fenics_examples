import fenics as fe
import numpy as np
from fenics import *
from mshr import *

#define the mesh
domain=Circle(Point(0,0),1)
mesh=generate_mesh(domain,64)
V=FunctionSpace(mesh,"P",1)

#define the pressure
alpha=4
R=0.6
beta=8
p=Expression("alpha*exp(-beta*(x[0]*x[0]+(x[1]-R)*(x[1]-R)))",degree=1,alpha=alpha,beta=beta,R=R)

#define dirchlet boundary
def boundary(x,on_boundary):
    return on_boundary
u_D=Constant(0)
bc=DirichletBC(V,u_D,boundary)

#define variational problem
u=TrialFunction(V)
v=TestFunction(V)
a=dot(grad(u),grad(v))*dx
L=v*p*dx

#solve the PDE problem
u=Function(V)
solve(a==L,u,bc)

#plot the solution
plot(u,title="Deflection")
#plot(mesh)

#store the solution in pvd file
vtkfile=File("membrane/solution.pvd")
vtkfile<<u

