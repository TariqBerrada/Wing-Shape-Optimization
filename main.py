import matplotlib.pyplot as plt
from MyStokes import *
from MyShapeOpt import *

parameters["reorder_dofs_serial"] = False

#Import the mesh.
omega = Mesh('stationary_circle.xml')
mf = MeshFunction('size_t', omega, 'stationary_circle_facet_region.xml')

#Get initial functions associated with mesh
[J, u, p, W, un] = StokesSolve(omega, mf)

#Define optimization parameters
n_iter = 20
rho = .003
epsilon = .001

#Optimize the mesh
mesh, J_list, u_optim, p_optim = optimize_mesh(omega, mf, n_iter, rho, epsilon)

print("Best loss obtained at iteration %d with a value of J=%s."%(J_list.index(min(J_list))+1, min(J_list)))

#Visualize the obtained solution
plt.figure()
plt.plot(J_list)
plt.title("Evolution of the loss function with each iteration.")
plt.xlabel("Iteration")
plt.ylabel("Loss: J(u)")


plt.figure()
plot(mesh)
plt.title("Optimized mesh.")

plt.figure()
plot(u_optim)
plt.title("u(x) at optimized mesh.")

plt.figure()
plot(p_optim)
plt.title("p(x) at optimized mesh.")

plt.show()


