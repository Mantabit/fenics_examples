from fenics import *
import sympy as sym
import fenics as fe

#define symbols for x and y and register them with sympy
x,y=sym.symbols("x[0],x[1]")

def q(u):
    return 1+u**2

#define sympy expressions for the boundary u and the inhomogenity f
u=1+x+2*y
f=-1*(sym.diff(q(u)*sym.diff(u,x),x)+sym.diff(q(u)*sym.diff(u,y),y))

f,u=sym.simplify(f),sym.simplify(u)

#compute the c code strings for the two expressuons f,u
f_ccode=sym.printing.ccode(f)
u_ccode=sym.printing.ccode(u)

print(f_ccode+"\n")
print(u_ccode+"\n")

#define the mesh and the trial an test function space
nx=ny=10
mesh=UnitSquareMesh(nx,ny)
V=FunctionSpace(mesh,"P",1)

#define the boundary and the boundary condition
def boundary(x,on_boundary):
    return on_boundary
u_D=Expression(u_ccode,degree=1)
bc=DirichletBC(V,u_D,boundary)

#define the variational problem
u=Function(V)
v=TestFunction(V)
f=Expression(f_ccode,degree=1)
F=f*v*dx+q(u)*dot(grad(v),grad(u))*dx

#solve the pde
solve(F==0,u,bc)

plot(u)