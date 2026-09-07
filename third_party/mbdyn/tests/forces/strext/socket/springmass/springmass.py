# $Header: /var/cvs/mbdyn/mbdyn/mbdyn-1.0/tests/forces/strext/socket/springmass/springmass.py,v 1.9 2017/01/12 15:03:07 masarati Exp $
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

import sys;
# set to path of MBDyn support for python communication
sys.path.append('/usr/local/mbdyn/libexec/mbpy');

import os;
import tempfile;
tmpdir = tempfile.mkdtemp('', '.mbdyn_');
path = tmpdir + '/mbdyn.sock';

os.environ['MBSOCK'] = path;
os.system('./mbdyn.sh -f springmass -o output > output.txt 2>&1 &');

from mbc_py_interface import mbcNodal
from numpy import *

host = "";
port = 0;
timeout = -1;	# forever
verbose = 0;
data_and_next = 1;
refnode = 0;	# no reference node
nodes = 1;	# number of nodes
labels = 0;
rot = 0x100;	# orientation vector
# rot = 0x200;	# orientation matrix
# rot = 0x400;	# Euler 123
accels = 1;
nodal = mbcNodal(path, host, port, timeout, verbose, data_and_next, refnode, nodes, labels, rot, accels);

nodal.negotiate();

nodal.n_f[0] = 0.;
nodal.n_f[1] = 0.;
nodal.n_f[2] = 0.;
nodal.n_m[0] = 0.;
nodal.n_m[1] = 0.;
nodal.n_m[2] = 0.;

# data
m = 1.;
omega = 10*2*pi;
xi = 0.1;
k = m*omega*omega;
d = 2*xi*m*omega;

i = 0;
maxiter = 5;
while 1:
	if (nodal.recv()):
		break;

	nodal.n_f[0] = -k*nodal.n_x[0] - d*nodal.n_xp[0];
	nodal.n_f[1] = -k*nodal.n_x[1] - d*nodal.n_xp[1];
	nodal.n_f[2] = -k*nodal.n_x[2] - d*nodal.n_xp[2];

	print 'i=', i, ' x=', nodal.n_x[2], ' xp=', nodal.n_xp[2], 'xpp=', nodal.n_xpp[2];

	# iterate until convergence
	i = i + 1;
	if (i == maxiter):
		i = 0;

	if (nodal.send(i == 0)):
		break;

nodal.destroy();
os.rmdir(tmpdir);

