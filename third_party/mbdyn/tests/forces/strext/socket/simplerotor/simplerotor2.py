from mbc_py_interface import mbcNodal
from numpy import *

# note: need mbc_py_interface.py and _mbc_py.so in PYTHONPATH; right now
# mbc_py_interface.py is in $(MBDYN_SRC)/libraries/libmbc
# _mbc_py.so is installed in $(MBDYN_libexecdir)/mbpy

path = "/tmp/mbdyn1.sock";
host = "";
port = 0;
timeout = -1;
verbose = 0;
data_and_next = 1;
refnode = 1;
nodes = 20;
labels = 0;
rot = 0x200;
accels = 0;
nodal1 = mbcNodal(path, host, port, timeout, verbose, data_and_next, refnode, nodes, labels, rot, accels);

path = "/tmp/mbdyn2.sock";
host = "";
port = 0;
timeout = -1;
verbose = 0;
data_and_next = 1;
refnode = 1;
nodes = 20;
labels = 0;
rot = 0x200;
accels = 0;
nodal2 = mbcNodal(path, host, port, timeout, verbose, data_and_next, refnode, nodes, labels, rot, accels);

nodal1.negotiate();
nodal2.negotiate();

fz = 1.e3;
nodal1.r_f[0] = 0.;
nodal1.r_f[1] = 0.;
nodal1.r_f[2] = 0.;
nodal1.r_m[0] = 0.;
nodal1.r_m[1] = 0.;
nodal1.r_m[2] = 0.;
nodal2.r_f[0] = 0.;
nodal2.r_f[1] = 0.;
nodal2.r_f[2] = 0.;
nodal2.r_m[0] = 0.;
nodal2.r_m[1] = 0.;
nodal2.r_m[2] = 0.;
for n in range(nodes):
	nodal1.n_f[3*n] = 0.;
	nodal1.n_f[3*n + 1] = 0.;
	nodal1.n_f[3*n + 2] = fz;
	nodal1.n_m[3*n] = 0.;
	nodal1.n_m[3*n + 1] = 0.;
	nodal1.n_m[3*n + 2] = 0.;

	nodal2.r_f[2] += nodal2.n_f[3*n + 2];

	nodal2.n_f[3*n] = 0.;
	nodal2.n_f[3*n + 1] = 0.;
	nodal2.n_f[3*n + 2] = fz;
	nodal2.n_m[3*n] = 0.;
	nodal2.n_m[3*n + 1] = 0.;
	nodal2.n_m[3*n + 2] = 0.;

	nodal2.r_f[2] += nodal2.n_f[3*n + 2];


while 1:
	if (nodal1.recv()):
		break;
	if (nodal2.recv()):
		break;

	# print 'nodal1: ', nodal1.r_omega[2];
	# print 'nodal2: ', nodal2.r_omega[2];
	print 'nodal1: ', nodal1.n_x[3*nodes-2];
	print 'nodal2: ', nodal2.n_x[3*nodes-2];

	if (nodal1.send(1)):
		break;
	if (nodal2.send(1)):
		break;

nodal1.destroy();
nodal2.destroy();

