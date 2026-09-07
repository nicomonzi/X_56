from __future__ import print_function
from mbc_py_interface import mbcNodal
from numpy import *

# note: need mbc_py_interface.py and _mbc_py.so in PYTHONPATH; right now
# mbc_py_interface.py is in $(MBDYN_SRC)/libraries/libmbc
# _mbc_py.so is installed in $(MBDYN_libexecdir)/mbpy

path = "/tmp/mbdyn.sock"
host = ""
port = 0
timeout = -1
verbose = 1
data_and_next = 1
refnode = 1
nodes = 40
labels = 0
rot = 0x100	# orientation vector
# rot = 0x200	# orientation matrix
# rot = 0x400	# Euler 123
accels = 0
nodal = mbcNodal(path, host, port, timeout, verbose, data_and_next, refnode, nodes, labels, rot, accels)

nodal.negotiate()

fz = 1.e3
nodal.r_f[0] = 0.
nodal.r_f[1] = 0.
nodal.r_f[2] = 0.
nodal.r_m[0] = 0.
nodal.r_m[1] = 0.
nodal.r_m[2] = 0.
for n in range(nodes):
	nodal.n_f[3*n] = 0.
	nodal.n_f[3*n + 1] = 0.
	nodal.n_f[3*n + 2] = fz
	nodal.n_m[3*n] = 0.
	nodal.n_m[3*n + 1] = 0.
	nodal.n_m[3*n + 2] = 0.

	nodal.r_f[2] += nodal.n_f[3*n + 2]
	

while 1:
	if (nodal.recv()):
		break;

	print('x = {}'.format(nodal.r_x));
	if (rot == 0x100):
		print('theta = {}'.format(nodal.r_theta))
	if (rot == 0x200):
		print('r = {} {} {}'.format(nodal.r_r[0], nodal.r_r[3], nodal.r_r[6]))
		print('{} {} {}'.format(nodal.r_r[1], nodal.r_r[4], nodal.r_r[7]))
		print('{} {} {}'.format(nodal.r_r[2], nodal.r_r[5], nodal.r_r[8]))
	if (rot == 0x400):
		print('e123 = {}'.format(nodal.r_euler_123))
	print('xp = {}'.format(nodal.r_xp))
	print('omega ='.format(nodal.r_omega));

	print('blade1 x = {}'.format(nodal.n_x[range(int(3*(nodes/2 - 1)), int(3*nodes/2))]))
	print('\t\ttheta = {}'.format(nodal.n_theta[range(int(3*(nodes/2 - 1)), int(3*nodes/2))]))
	print('blade2 x = {}'.format(nodal.n_x[range(int(3*(nodes - 1)), int(3*nodes))]))
	print('\t\ttheta = {}'.format(nodal.n_theta[range(int(3*(nodes - 1)), int(3*nodes))]))

	if (nodal.send(1)):
		break;

nodal.destroy();

