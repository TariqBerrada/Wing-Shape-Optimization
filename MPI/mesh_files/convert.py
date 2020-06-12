from fenics import *

for k in range(1, 5):
	omega = Mesh("mesh_{}.xml".format(k))
	mf = MeshFunction('size_t', omega, 'mesh_{}_facet_region.xml'.format(k))
	hdf = HDF5File(omega.mpi_comm(), "mesh_{}.h5".format(k), "w")
	hdf.write(omega, "/mesh")
	hdf.write(mf, "/regions")