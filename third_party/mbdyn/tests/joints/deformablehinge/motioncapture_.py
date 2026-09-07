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
# exec "mbdyn body"
# exec "python spring.py"
from __future__ import print_function
import socket
import struct
import sys;
from numpy import *
import os;
import tempfile;
import time;

tmpdir = tempfile.mkdtemp('', '.mbdyn_');
path_in = tmpdir + '/mbdyn_in.sock';
path_out = tmpdir + '/mbdyn_out.sock';

os.environ['MBSKIN'] = path_in;
os.environ['MBSKOU'] = path_out;
os.system('./mbdyn.sh -f motioncapture -o output -N 1 > output.txt 2>&1 &');

# wait for MBDyn to start and create sockets...
time.sleep(2);

# create input socket (4 double: Fx, Fy of 2 sensors)
# NOTE: what's out for MBDyn is in for this script
print('input socket: %s' % path_out);
s_in = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM);
s_in.connect(path_out);
s_in_bufsize = 4*8;

# create output socket (4 double: eps_x, eps_y of 2 sensors)
# NOTE: what's in for MBDyn is out for this script
print('output socket: %s' % path_in);
s_out = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM);
s_out.connect(path_in);
s_out_bufsize = 4*8;

# stiffness
K = 1.e3;
alpha = 0.5;

# initial value
eps_x_1 = 0.;
eps_y_1 = 0.;
eps_x_2 = 0.;
eps_y_2 = 0.;

f_x_1 = 0.;
f_y_1 = 0.;
f_x_2 = 0.;
f_y_2 = 0.;

# step counter
i = 0;
while 1:
	# read motion (x, x_prime) from body's output element
	buf_in = bytearray(s_in_bufsize);
	rc = s_in.recv_into(buf_in, s_in_bufsize);
	if (rc != s_in_bufsize):
		print('received %d bytes; closing...' % rc);
		break;

	f_x_1, f_y_1, f_x_2, f_y_2 = struct.unpack('dddd', buf_in);

	# compute eps for next step
	eps_x_1 = eps_x_1 + alpha*f_x_1/K;
	eps_y_1 = eps_y_1 + alpha*f_y_1/K;
	eps_x_2 = eps_x_2 + alpha*f_x_2/K;
	eps_y_2 = eps_y_2 + alpha*f_y_2/K;

	print('i= %d f= %e %e %e %e  eps= %e %e %e %e' % (i, f_x_1, f_y_1, f_x_2, f_y_2, eps_x_1, eps_y_1, eps_x_2, eps_y_2));

	# send force to body's file driver
	buf_out = bytearray(struct.pack('dddd', eps_x_1, eps_y_1, eps_x_2, eps_y_2));
	#rc = s_out.send(buf_out, s_out_bufsize)
	rc = s_out.send(buf_out);
	if (rc != s_out_bufsize):
		print('sent %d bytes; closing...' % rc);
		break;

	i = i + 1;

s_out.close();
s_in.close();

os.rmdir(tmpdir);
