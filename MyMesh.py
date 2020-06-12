from fenics import *
from dolfin import *
import numpy as np
import matplotlib.pyplot as plt

def vertex_normal(mesh, boundary_parts, i):
    """Calculates the normal in a given vertex of the mesh."""
    ver = Vertex(mesh, i)
    n = Point(0.0, 0.0)
    div = 0.0
    for fac in entities(ver, mesh.geometry().dim()-1):
        f = Facet(mesh, fac.index())
        if f.exterior()==True:
            n+=f.normal()
            div+=1
    n/=div
    return n

def grad2n(nx, ny, u, un):
    """returns the functions representing the x and y component of gard(u).grad(u).n"""
    Vn = u.function_space()
    mesh = Vn.mesh()
    Wn = FunctionSpace(mesh, 'CG', 1)

    grad_nu_x =inner(un, nx)
    grad_nu_y = inner(un, ny)

    grad_nu_x = project(grad_nu_x, Wn)
    grad_nu_y = project(grad_nu_y, Wn)

    gnux = Function(Wn, grad_nu_x.vector())
    gnuy = Function(Wn, grad_nu_y.vector())
    return gnux, gnuy

def compute_normal_field(omega, mf):
    """Computes a vector field whose elements are the normals at each point."""
    #Boundary Mesh
    gamma = BoundaryMesh(omega, "exterior")
    mapa = gamma.entity_map(0)
    
    V = FunctionSpace(gamma, "CG", 1)
    N = VectorFunctionSpace(gamma, "CG", 1)
    normal_field = Function(N)
    
    normal_x = Function(V)
    normal_y = Function(V)
    
    for ver in entities(gamma, 0):
        i = mapa[ver.index()]

        point = vertex_normal(omega, mf, i)
        normal_x.vector()[ver.index()] = -1.0*point[0]
        normal_y.vector()[ver.index()] = -1.0*point[1]
    
    (x, y) = TestFunctions(N)
    solve(inner(normal_field[0], x)*dx + inner(normal_field[1], y)*dx - inner(normal_x, x)*dx -inner(normal_y, y)*dx == 0, normal_field)    
	
    #Dirichlet Boundary must be initialized by Function in Vector Space over Omega (only on Gamma fails)
    V_vec = VectorFunctionSpace(omega, "CG", 1)
    normal_fieldV = Function(V_vec)
    for ver in entities(gamma, 0):
        i = mapa[ver.index()]
        normal_fieldV.vector()[i] = normal_field.vector()[ver.index()]
        normal_fieldV.vector()[i+omega.num_vertices()] = normal_field.vector()[ver.index()+gamma.num_vertices()]
    
    deform = TrialFunction(V_vec)
    v = TestFunction(V_vec)
    a = 0.01*inner(nabla_grad(deform), nabla_grad(v))*dx + 0.01*inner(deform,v)*ds(4)
    L = inner(Constant((0.0,0.0)),v)*dx
    bc1 = DirichletBC(V_vec, normal_fieldV, mf, 15)     #Obstacle
    bc2 = DirichletBC(V_vec, Constant((0,0)), mf, 13)    #outflow
    bc3 = DirichletBC(V_vec, Constant((0,0)), mf,14)    #noslip
    bc4 = DirichletBC(V_vec, Constant((0,0)), mf, 12)    #inflow
    bc = [bc1, bc2, bc3, bc4]
    deform = Function(V_vec)
    solve(a==L, deform, bcs=bc)
    return deform

def mesh_surface(mesh):
	"""Calculates surface of a mesh."""
	Wu = FunctionSpace(mesh, 'CG', 1)
	e = Expression('1.0', degree=1)
	unity = interpolate(e, Wu)
	surface = unity*dx
	surface = assemble(surface)
	return surface