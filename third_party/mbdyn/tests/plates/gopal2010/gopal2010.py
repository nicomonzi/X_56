#! /usr/bin/python
#
# $Header: /var/cvs/mbdyn/mbdyn/mbdyn-1.0/tests/plates/gopal2010/gopal2010.py,v 1.2 2014/12/13 17:15:55 masarati Exp $
#
# MBDyn (C) is a multibody analysis code. 
# http://www.mbdyn.org
# 
# Copyright (C) 1996-2010
# 
# Pierangelo Masarati	<masarati@aero.polimi.it>
# Paolo Mantegazza	<mantegazza@aero.polimi.it>
# 
# Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
# via La Masa, 34 - 20156 Milano, Italy
# http://www.aero.polimi.it
# 
# Changing this copyright notice is forbidden.
# 
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation (version 2 of the License).
# 
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
# 

import sys;
import os;
import tempfile;

sys.path.append(os.path.abspath('/usr/local/mbdyn/libexec/mbpy/'));

from mbc_py_interface import mbcNodal
from numpy import *

# note: need mbc_py_interface.py, mbc_py.py and _mbc_py.so in PYTHONPATH
# right now they are all installed in $(MBDYN_libexecdir)/mbpy
# execute after setting PYTHONPATH=$(MBDYN_libexecdir)/mbpy

tmpdir = tempfile.mkdtemp('', '.mbdyn_');
path = tmpdir + '/mbdyn.sock';
os.environ['MBSOCK'] = path;

#inputfile = 'gopal2010_ext';
inputfile = 'gopal2010m_ext';

mbdyn_cmd = 'mbdyn';
os.system(mbdyn_cmd + ' -f ' + inputfile + ' -o output > output.txt 2>&1 &');

host = "";
port = 0;
timeout = -1;
verbose = 0;
data_and_next = 1;
refnode = 1;
nodes = 34595 # set to the actual count of mapped points
labels = 0;
# rot = 0x1000;	# reference node orientation vector; no node orientations
rot = 0x2000;	# reference node orientation matrix; no node orientations
# rot = 0x4000;	# reference node Euler 123; no node orientations
accels = 0;
structure = mbcNodal(path, host, port, timeout, verbose, data_and_next, refnode, nodes, labels, rot, accels);

# note: could fail!
structure.negotiate();

# structure.r_*: reference node stuff (if configured)
# structure.n_*: node stuff

# input:
# structure.r_x: reference node position
# structure.r_theta: reference node orientation vector (if rot == 0x1000)
# structure.r_r: reference node orientation matrix (if rot == 0x2000)
# structure.r_euler_123: reference node Euler angles 123 (if rot == 0x4000)
# structure.r_xp: reference node velocity
# structure.r_omega: reference node angular velocity
# structure.r_xpp: reference node acceleration (if accels != 0)
# structure.r_omegap: reference node angular acceleration (if accels != 0)

# output:
# structure.r_f: reference node force
# structure.r_m: reference node moment

# input:
# structure.n_x: nodes position
# structure.n_xp: nodes velocity
# structure.n_xpp: nodes acceleration (if accels != 0)

# output:
# structure.n_f: nodes force

# examples:
# position of reference node: structure.r_x[0], structure.r_x[1], structure.r_x[2]
# position of node 7: structure.n_x[7*3-3], structure.n_x[7*3-2], structure.n_x[7*3-1]

# more than one communicator can be used in the same analysis

structure.r_f[0] = 0.;
structure.r_f[1] = 0.;
structure.r_f[2] = 0.;
structure.r_m[0] = 0.;
structure.r_m[1] = 0.;
structure.r_m[2] = 0.;

#xmprev = zeros((3*nodes, 1));
#xmcurr = zeros((3*nodes, 1));
#xm = zeros((3*nodes, 1));

# damp = -2e-8;
damp = -2e-7;
# damp = -2e-6;
# damp = -2e-5;
cnt = 0;
iter_cnt = 0;
while 1:
	if (structure.recv()):
		break;

	x0 = reshape(structure.r_x[0:3], (3, 1));
	R0 = reshape(structure.r_r, (3, 3), 'F');
	x0p = dot(R0.transpose(), structure.r_xp);
	omega0 = dot(R0.transpose(), structure.r_omega);

	#for n in range(nodes):
	#	xmcurr[3*n:3*(n+1)] = dot(R0, reshape(structure.n_x[3*n:3*(n+1)], (3, 1)));
	#if (cnt > 0):
	#	xm = reshape(xmcurr, (3*nodes, 1)) - reshape(xmprev, (3*nodes, 1));
	#xmprev = xmcurr.copy();

	print '#'
	#xm = ravel(xm);
	#for n in range(nodes):
	#	print xm[3*n], xm[3*n+1], xm[3*n+2]

	for n in range(nodes):
		structure.n_f[3*n:3*(n+1)] = damp*(x0p + cross(omega0, structure.n_x[3*n:3*(n+1)]) + structure.n_xp[3*n:3*(n+1)]);

	if iter_cnt == 5:
		print '# cnt', cnt
		for n in range(nodes):
			x = ravel(x0 + dot(R0, reshape(structure.n_x[3*n:3*(n+1)], (3, 1))));
			print x[0], x[1], x[2]
			# print structure.n_x[3*n], structure.n_x[3*n + 1], structure.n_x[3*n + 2];
	
	# send(1) at convergence; send(0) for another iteration
	if (structure.send(iter_cnt == 5)):
		break;

	if iter_cnt == 5:
		iter_cnt = 0;
		cnt = cnt + 1;

	iter_cnt = iter_cnt + 1;

structure.destroy();
os.rmdir(tmpdir);

