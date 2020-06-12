# Script python expérimental qui a pour but de générer automatiquement des sous-maillages du maillage principal

from os.path import join
from fenics import *
from MyStokes import *
import matplotlib.pyplot as plt


# Import the mesh.
omega = Mesh('stationary_circle.xml')
mf = MeshFunction('size_t', omega, 'stationary_circle_facet_region.xml')


def createSubMeshes(nb_thread, omega, overlap=0.2):
	"""Return a list of nb_thread submeshes."""

	tol = 1E-14
	res = []

	X = []
	for ver in vertices(omega):
		X.append(ver.point().x())

	X.sort()
	overlap = int((overlap * len(X)) // nb_thread)

	domain = MeshFunction("size_t", omega, omega.topology().dim())

	for k in range(nb_thread):
		min_x = max(0, k * len(X)//nb_thread - overlap)
		max_x = min(len(X) - 1, (k + 1) * len(X)//nb_thread + overlap)  # Il peut y avoir des problèmes si nb_thread * len(X)//nb_thread + overlap n'atteint pas len(X) - 1

		min_x = X[min_x]
		max_x = X[max_x]

		desc = "min_x - tol <= x[0] && x[0] <=  max_x + tol"

		subdomain = CompiledSubDomain(desc, tol=tol, min_x=min_x, max_x=max_x)

		domain.set_all(0)  # Remise à 0 des compteurs
		subdomain.mark(domain, 1)
		submesh = SubMesh(omega, domain, 1)

		res.append(submesh)

	return res


def globalIndex(submesh, i):
	return submesh.data().array('parent_vertex_indices', 0)[i]


def determineBorder(submeshes):
	"""
	Return a list of list
	Each list contain the list of the boundaries of the others submeshes which are in the submesh
	"""
	n = len(submeshes)
	index_border = [[] for _ in range(n)]
	for k in range(n):
		b = BoundaryMesh(submeshes[k], "exterior")
		mapb = b.entity_map(0)
		for ver in entities(b, 0):
			i = mapb[ver.index()]
			i = globalIndex(submeshes[k], i)
			index_border[k].append(i)

	index_all = [[] for _ in range(n)]
	for k in range(n):
		m = submeshes[k]
		for ver in vertices(m):
			i = globalIndex(m, ver.index())
			index_all[k].append(i)

	res = [[] for _ in range(n)]
	for s in range(n):
		for t in range(n):
			for i in index_border[s]:
				if i in index_all[t]:
					res[t][s].append(i)
	return res






nb_thread = 4
meshes = createSubMeshes(nb_thread, omega)

# print(meshes[0].data().array('parent_vertex_indices', 0)[0])
# print(determineBorder(meshes))


# mesh0 = meshes[0]
# mf0 = MeshFunction('size_t', mesh0, 'stationary_circle_facet_region.xml')
# [J, u, p, W, un] = StokesSolve(mesh0, mf0)

for t in range(nb_thread):
	plt.clf()
	plot(omega)
	plot(meshes[t], color="blue")
	plt.savefig(join("submeshes", "mesh{}.png".format(t)))
	plt.clf()
	plot(meshes[t])
	plt.savefig(join("submeshes", "detail", "mesh{}.png".format(t)))
