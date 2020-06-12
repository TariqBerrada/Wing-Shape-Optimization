# import matplotlib.pyplot as plt
from MyStokes import *
from MyShapeOpt import *
from fenics import *
import time
from mpi4py import MPI


def collides(mesh, point):
	return mesh.bounding_box_tree().compute_first_entity_collision(point) < 2 *omega.num_vertices()


def update_solution(mesh, mf, u, u2):
	#[_, u, _, _, _] = StokesSolve(mesh, mf)
	V = FunctionSpace(mesh, "CG", 1)
	N = VectorFunctionSpace(mesh, "CG", 1)

	u_final = Function(N)
	u_dW_x = Function(V)
	u_dW_y = Function(V)
	for vertex in vertices(mesh):
		point = vertex.point()
		point_o = u(vertex.point())

		Vu = u2.function_space()
		other_mesh = Vu.mesh()

		if collides(other_mesh, point):
			u_dW_x.vector()[vertex.index()] = u2(point)[0]
			u_dW_y.vector()[vertex.index()] = u2(point)[1]
		else:
			u_dW_x.vector()[vertex.index()] = point_o[0]
			u_dW_y.vector()[vertex.index()] = point_o[1]

	(x, y) = TestFunctions(N)
	solve(inner(u_final[0], x)*dx + inner(u_final[1], y)*dx - inner(u_dW_x, x)*dx - inner(u_dW_y, y)*dx == 0, u_final)

	return u_final


def optimization(u_list, mesh_list, mf_list, epochs):
	"""Old sequential function"""
	u_new_list = u_list
	for epoch in range(epochs):
		print("Epoch %d out of %d"%(epoch+1, epochs))
		for i in range(len(u_list)):
			for j in range(len(u_list)):
				if j != i:
					u_new_list[i] = update_solution(mesh_list[i], mf_list[i],u_new_list[i], u_list[j])
	return u_new_list


def parallel_optimization(mesh, mf, epochs):
	"""Same that optimization, but in parallel"""
	# Get initial functions associated with mesh
	[_, u, _, _, _] = StokesSolve(mesh, mf)
	u_list = echange_u(u)

	for epoch in range(epochs):
		if Me == 0:
			print("Epoch %d out of %d"%(epoch+1, epochs))
		for other_u in u_list:
			u = update_solution(mesh, mf, u, other_u)
		u_list = echange_u(u)
	return u


def echange_u(u):
	ur = u
	u_list = []

	if Me % 2 == 0:
		comm.Ssend(u, (Me+1)%NbP)

		comm.Recv(ur, source=(Me-1)%NbP)
		u_list.append(up)

		comm.Ssend(u, (Me-1)%NbP)

		comm.Recv(ur, source=(Me+1)%NbP)
		u_list.append(ur)

		comm.Ssend(u, (Me+2)%NbP)

		comm.Recv(ur, source=(Me-2)%NbP)
		u_list.append(ur)

	else:
		comm.Recv(ur, source=(Me-1)%NbP)
		u_list.append(ur)

		comm.Ssend(u, (Me+1)%NbP)

		comm.Recv(ur, source=(Me+1)%NbP)
		u_list.append(ur)

		comm.Ssend(u, (Me-1)%NbP)

		comm.Recv(ur, source=(Me-2)%NbP)
		u_list.append(ur)

		comm.Ssend(u, (Me+2)%NbP)

	return u_list




########################################################################
# Main function
########################################################################

# MPI information extraction
comm = MPI.COMM_WORLD
NbP = comm.Get_size()
Me  = comm.Get_rank()


epochs = 3

parameters["reorder_dofs_serial"] = False
omega = Mesh('./mesh_files/mesh.xml')


if NbP != 4:
	print(NbP, "process détécté. Il en faut que 4.")
else:
	tpar = time.time()

	meshname = "/mnt/d/CentraleSupelec/Calcul haute performance/ONERA/project/parallel/mesh_files/mesh_{}.h5".format(Me+1)
	omegaMe = Mesh()
	hdf = HDF5File(omegaMe.mpi_comm(),  meshname, "r")
	hdf.read(omegaMe, "/mesh", False)
	mfMe = MeshFunction('size_t', omegaMe, omega.topology().dim())
	hdf.read(mfMe, "/region", False)

	# omegaMe = Mesh("./mesh_files/mesh_{}.xml".format(Me+1))
	# mfMe = MeshFunction('size_t', omegaMe, './mesh_files/mesh_{}_facet_region.xml'.format(Me+1))

	u = parallel_optimization(omegaMe, mfMe, epochs)
	# u is not a list! Each process have is own u. The former u_new can be obtain by gathering all the u over the different process
	print("\nTemps d'execution en parallèle :", time.time()-tpar)


if Me == 0:
	## ===========================================
	# Former varable used in optimization
	tsec = time.time()
	omega1 = Mesh('./mesh_files/mesh_1.xml')
	omega2 = Mesh('./mesh_files/mesh_2.xml')
	omega3 = Mesh('./mesh_files/mesh_3.xml')
	omega4 = Mesh('./mesh_files/mesh_4.xml')
	mf1 = MeshFunction('size_t', omega1, './mesh_files/mesh_1_facet_region.xml')
	mf2 = MeshFunction('size_t', omega2, './mesh_files/mesh_2_facet_region.xml')
	mf3 = MeshFunction('size_t', omega3, './mesh_files/mesh_3_facet_region.xml')
	mf4 = MeshFunction('size_t', omega4, './mesh_files/mesh_4_facet_region.xml')

	omega = Mesh('./mesh_files/mesh.xml')

	parameters["reorder_dofs_serial"] = False

	# Get initial functions associated with mesh
	[J1, u1, p1, W1, un1] = StokesSolve(omega1, mf1)
	[J2, u2, p2, W2, un2] = StokesSolve(omega2, mf2)
	[J3, u3, p3, W3, un3] = StokesSolve(omega3, mf3)
	[J4, u4, p4, W4, un4] = StokesSolve(omega4, mf4)

	u_list = [u1, u2, u3, u4]
	mesh_list = [omega1, omega2, omega3, omega4]
	mf_list = [mf1, mf2, mf3, mf4]

	## ===========================================
		

	u_new = optimization(u_list, mesh_list, mf_list, epochs)

	print("\nTemps d'execution en sequentiel :", time.time()-tsec)
