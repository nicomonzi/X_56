#! /usr/bin/python
#
# $Header: /var/cvs/mbdyn/mbdyn/mbdyn-1.0/tests/plates/agard/agard_mapping.py,v 1.8 2017/01/12 15:06:33 masarati Exp $
#
# MBDyn (C) is a multibody analysis code. 
# http://www.mbdyn.org
# 
# Copyright (C) 1996-2017
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

from __future__ import print_function
from mbc_py_interface import mbcNodal
from numpy import *

# note: need mbc_py_interface.py and _mbc_py.so in PYTHONPATH; right now
# mbc_py_interface.py is in $(MBDYN_SRC)/libraries/libmbc
# _mbc_py.so is installed in $(MBDYN_libexecdir)/mbpy

path = "/tmp/mbdyn_agard.sock"
host = ""
port = 0
timeout = -1
verbose = 0
data_and_next = 1
refnode = 1
nodes = 1050	# set to the actual count of mapped points
# nodes = 12025	# set to the actual count of mapped points
labels = 0
# rot = 0x1000	# reference node orientation vector; no node orientations
rot = 0x2000	# reference node orientation matrix; no node orientations
# rot = 0x4000	# reference node Euler 123; no node orientations
accels = 0
nodal = mbcNodal(path, host, port, timeout, verbose, data_and_next, refnode, nodes, labels, rot, accels)

nodal.negotiate()

fz = 1.e-1
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

	nodal.r_f[2] += nodal.n_f[3*n + 2]
	
test = 1e9
tol = 1e-9
old_x = 0.
while 1:
	if (nodal.recv()):
		break

	print('x121=\t{}'.format(nodal.n_x[3*(nodes - 1) + 0], nodal.n_x[3*(nodes - 1) + 1], nodal.n_x[3*(nodes - 1) + 2]))

	# simple convergence test based on z component of last node
	new_x = nodal.n_x[3*(nodes - 1) + 2]
	test = abs(new_x - old_x)
	old_x = new_x
	res = 0
	if (test < tol):
		res = 1

	print('test={} res={} '.format(test,res));

	# send(1) at convergence; send(0) for another iteration
	if (nodal.send(res)):
		break;

nodal.destroy();

