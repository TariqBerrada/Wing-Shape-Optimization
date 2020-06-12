from MyStokes import *
from MyShapeOpt import *
from fenics import *
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


epochs = 3

parameters["reorder_dofs_serial"] = False
omega = Mesh('./mesh_files/mesh.xml')


omega1 = Mesh('./mesh_files/mesh_1.xml')
omega2 = Mesh('./mesh_files/mesh_2.xml')
omega3 = Mesh('./mesh_files/mesh_3.xml')
omega4 = Mesh('./mesh_files/mesh_4.xml')
mf1 = MeshFunction('size_t', omega1, './mesh_files/mesh_1_facet_region.xml')
mf2 = MeshFunction('size_t', omega2, './mesh_files/mesh_2_facet_region.xml')
mf3 = MeshFunction('size_t', omega3, './mesh_files/mesh_3_facet_region.xml')
mf4 = MeshFunction('size_t', omega4, './mesh_files/mesh_4_facet_region.xml')


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
