from fenics import *
from dolfin import *
import numpy as np
import matplotlib.pyplot as plt
from MyMesh import *
from MyStokes import *

def gradient_params(J_list, rho):
    """Use this argument to have a variable parameter in gradient descnet"""
    grad_J = np.abs(J_list[-1] - J_list[-2])
    rho1 = rho
    if grad_J>.02:
    	rho1 = 4*rho
    if grad_J < .02:
    	rho1 = rho/2
    return rho1


def update_mesh(omega, u, gradx, grady, rho, initial_surface, nx, ny, epsilon):
    """moves the vertices of the mesh according to the direction of the gradient."""

    omega.coordinates()[:, 0] -= rho*np.asarray(gradx.vector())
    omega.coordinates()[:, 1] -= rho*np.asarray(grady.vector())

    Wn = FunctionSpace(omega, 'CG', 1)

    dir_x = project(nx, Wn)
    dir_y = project(ny, Wn)

    while np.abs(mesh_surface(omega) - initial_surface) > .02:
        omega.coordinates()[:, 0] += epsilon*np.asarray(dir_x.vector())
        omega.coordinates()[:, 1] += epsilon*np.asarray(dir_y.vector())
        
    return omega

def optimize_mesh(mesh, mf, n_iter, rho, epsilon):
    """Does gradient descent in order to optimize the shape of the mesh. """
    omega = mesh
    initial_surface = mesh_surface(mesh)
    J_list = []
    filename = "./evolution/mesh/mesh_at_iter_0.png"
    plt.figure()
    plot(omega)
    plt.savefig(filename)
    plt.close()
    u_optim = None
    p_optim = None
    J_optim = np.inf
    mesh_optim = None

    for i in range(n_iter):
        print("Iteration " + str(i+1)+" out of " + str(n_iter) + ".")
        [J, u, p, W, un] = StokesSolve(omega, mf)
        J_list.append(J)

        if J < J_optim:
        	J_optim = J
        	u_optim = u
        	p_optim = p
        	mesh_optim = omega

        #Calculate the normal vectors
        n_field = compute_normal_field(omega, mf)
        nx, ny = n_field.sub(0), n_field.sub(1)

        #Calculate the gradient vectors
        gnux, gnuy = grad2n(nx, ny, u, un)

        #Update the mesh accordingly
        omega = update_mesh(omega, u, gnux, gnuy, rho, initial_surface, nx, ny, epsilon)

        #Saving the mesh frame
        filename = "./evolution/mesh/mesh_at_iter_"+str(i+1)+".png"
        plt.figure()
        plot(omega)
        plt.savefig(filename)
        plt.close()

        #Saving u
        filename_u = "./evolution/u/mesh_at_iter_"+str(i+1)+".png"
        plt.figure()
        plot(u)
        plt.savefig(filename_u)
        plt.close()

        #Saving p 
        filename_p = "./evolution/p/mesh_at_iter_"+str(i+1)+".png"
        plt.figure()
        plot(p)
        plt.savefig(filename_p)
        plt.close()        
    return mesh_optim, J_list, u_optim, p_optim

