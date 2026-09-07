from mbc_py_interface import mbcNodal
from numpy import *

# note: need mbc_py_interface.py and _mbc_py.so in PYTHONPATH; right now
# mbc_py_interface.py is in $(MBDYN_SRC)/libraries/libmbc
# _mbc_py.so is installed in $(MBDYN_libexecdir)/mbpy

path = "/tmp/mbdyn.sock";
host = "";
port = 0;
timeout = -1;
verbose = 0;
data_and_next = 1;
refnode = 1;
nodes = 40;
labels = 0;
rot = 0x000;	# none
# rot = 0x100;	# orientation vector
# rot = 0x200;	# orientation matrix
# rot = 0x400;	# Euler 123
accels = 0;
nodal = mbcNodal(path, host, port, timeout, verbose, data_and_next, refnode, nodes, labels, rot, accels);

nodal.negotiate();

fz = 1.e3;
nodal.r_f[0] = 0.;
nodal.r_f[1] = 0.;
nodal.r_f[2] = 0.;
#nodal.r_m[0] = 0.;
#nodal.r_m[1] = 0.;
#nodal.r_m[2] = 0.;
for n in range(nodes):
	nodal.n_f[3*n] = 0.;
	nodal.n_f[3*n + 1] = 0.;
	nodal.n_f[3*n + 2] = fz;
	#nodal.n_m[3*n] = 0.;
	#nodal.n_m[3*n + 1] = 0.;
	#nodal.n_m[3*n + 2] = 0.;

	nodal.r_f[2] += nodal.n_f[3*n + 2];
	

while 1:
	if (nodal.recv()):
		break;

	print 'x=    ', nodal.r_x;
	if (rot == 0x100):
		print 'theta=', nodal.r_theta;
	if (rot == 0x200):
		print 'r=    ', nodal.r_r[0], nodal.r_r[3], nodal.r_r[6];
		print '      ', nodal.r_r[1], nodal.r_r[4], nodal.r_r[7];
		print '      ', nodal.r_r[2], nodal.r_r[5], nodal.r_r[8];
	if (rot == 0x400):
		print 'e123= ', nodal.r_euler_123;
	print 'xp=   ', nodal.r_xp;
	#print 'omega=', nodal.r_omega;

	print 'blade1 x=    ', nodal.n_x[range(3*(nodes/2 - 1), 3*nodes/2)];
	#print '       theta=', nodal.n_theta[range(3*(nodes/2 - 1), 3*nodes/2)];
	print 'blade2 x=    ', nodal.n_x[range(3*(nodes - 1), 3*nodes)];
	#print '       theta=', nodal.n_theta[range(3*(nodes - 1), 3*nodes)];

	if (nodal.send(1)):
		break;

nodal.destroy();

