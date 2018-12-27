import fenics as fe
import numpy as np
from fenics import *

#define the function space by specifying the mesh and the order of approximation
mesh=fe.UnitSquareMesh(8,8)
V=fe.FunctionSpace(mesh,"P",1)

#define an expression which will define our Dirichlet boundary condition
u_D=fe.Expression("1+pow(x[0],2)+2*pow(x[1],2)",degree=2)

#define a boolean function which describes whether a point is on a Dirichlet boundary
def boundary(x,on_boundary):
    return on_boundary
bc=fe.DirichletBC(V,u_D,boundary)

#define variational problem
u=fe.TrialFunction(V)
v=fe.TestFunction(V)
f=fe.Constant(-6.0)
a=(fe.dot(fe.grad(u),fe.grad(v)))*fe.dx
L=f*v*fe.dx

#compute the solution
u=fe.Function(V)
fe.solve(a==L,u,bc)

#plot solution and mesh
fe.plot(u)
fe.plot(mesh)

#save the solution into a vtk file
vtkfile=fe.File("poisson/solution.pvd")
vtkfile<<u

#compute the error
error_l2=fe.errornorm(u_D,u,"L2")
print("The overall L2 error is given by: "+str(error_l2)+"\n")

#maximum error at vertices
vertex_values_u_D=u_D.compute_vertex_values(mesh)
vertex_values_u=u.compute_vertex_values(mesh)
max_vertex_error=np.max(np.abs(vertex_values_u_D-vertex_values_u))
print("Maximum vertex error given by: "+str(max_vertex_error)+"\n")