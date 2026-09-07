# $Header: /var/cvs/mbdyn/mbdyn/mbdyn-1.0/tests/forces/modalext/socket/socket.py,v 1.5 2017/01/12 15:02:26 masarati Exp $
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

import sys
# set to path of MBDyn support for python communication
sys.path.append('/usr/local/mbdyn/libexec/mbpy')

import os
import tempfile
tmpdir = tempfile.mkdtemp('', '.mbdyn_')
path = tmpdir + '/mbdyn.sock'
#path = '/tmp/mbdyn.sock'

os.environ['MBSOCK'] = path
os.system('./mbdyn.sh -f socket -o output > output.txt 2>&1 &')

from mbc_py_interface import mbcModal
from numpy import *

host = ""
port = 0
timeout = -1	# forever
verbose = 0
data_and_next = 1
refnode = 1	# modal node is present, rigid body motion is dealt with
modes = 4	# number of modes
modal = mbcModal(path, host, port, timeout, verbose, data_and_next, refnode, modes)

modal.negotiate()

# initialize rigid and modal forces
modal.r_f[0] = 0.
modal.r_f[1] = 0.
modal.r_f[2] = 0.
modal.r_m[0] = 0.
modal.r_m[1] = 0.
modal.r_m[2] = 0.

modal.m_p[0] = 0.
modal.m_p[1] = 0.
modal.m_p[2] = 0.
modal.m_p[3] = 0.

# data (simple mass/spring/damper for all modes)
m = 1.
omega = 10*2*pi
xi = 0.1
k = m*omega*omega
d = 2*xi*m*omega

step = 0
i = 0
maxiter = 5
while 1:
	if (modal.recv()):
		break

	if (step == 0):
		dp = 1. # perturbation
	else:
		dp = 0.

	modal.r_f[0] = 0.*modal.r_x[0]
	modal.r_f[1] = 0.*modal.r_x[1]
	modal.r_f[2] = 0.*modal.r_x[2]

	modal.m_p[0] = dp - k*modal.m_q[0] - d*modal.m_qp[0]
	modal.m_p[1] = dp - k*modal.m_q[1] - d*modal.m_qp[1]
	modal.m_p[2] = dp - k*modal.m_q[2] - d*modal.m_qp[2]
	modal.m_p[3] = dp - k*modal.m_q[3] - d*modal.m_qp[3]

	# iterate until convergence
	i = i + 1;
	if (i == maxiter):
		i = 0
		step = step + 1

	if (modal.send(i == 0)):
		break

modal.destroy()
os.rmdir(tmpdir)

