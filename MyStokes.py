# -*- coding: utf-8 -*-
from fenics import *
from dolfin import *
import numpy as np
import matplotlib.pyplot as plt

def step(u):
	"""Function used in Stokes Solve"""
	return inner(grad(u), grad(u))

def StokesSolve(omega, mf):
	"""Solves the Stokes physical problem on a given mesh."""

	V = VectorElement( "CG", omega.ufl_cell(), 2)
	Q = FiniteElement("CG", omega.ufl_cell(), 1)
	W = FunctionSpace(omega, MixedElement(V,Q))

	(u, p) = TrialFunctions(W)
	(v, q) = TestFunctions(W)

	f1 = Constant((0.0, 0.0))

	a = inner(nabla_grad(u),nabla_grad(v))*dx + p*div(v)*dx - div(u)*q*dx
	L = inner(f1,v)*dx

	#TP4 boundary Conditions for stationary circle.
	bc1 = DirichletBC(W.sub(0), (0.0, 0.0), mf, 14)
	bc2 = DirichletBC(W.sub(0), (0.0, 0.0), mf, 15)    #Obstacle
	bc3 = DirichletBC(W.sub(0), (1.0, 0.0), mf, 12)    #Inflow

	bc = [bc1, bc2, bc3]

	#Calculate u and p.
	w = Function(W)
	solve (a == L, w, bc)
	(u,p) = w.split(True)
    
	#Calculate loss.
	J = 0.5*inner(grad(u),grad(u))*dx
	J = assemble(J)

	Vn = u.function_space()
	mesh = Vn.mesh()
	
	#Calculate grad(u).grad(u).
	Wn = FunctionSpace(mesh, 'CG', 1)
	un = step(u)
	un = project(un, Wn)
	
	return J, u, p,  W, un


