# MBDyn (C) is a multibody analysis code.
# http://www.mbdyn.org
#
# Copyright (C) 1996-2024
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
# Cosimulation output content module test

import socket
import struct
import os;
import tempfile;
import time;

tmpdir = tempfile.mkdtemp('', '.mbdyn_');
path_out = tmpdir + '/mbdyn_out.sock';

os.environ['MBSKOU'] = path_out;

try:
        gtest_mbdyn_args = os.environ['GTEST_MBDYN_ARGS'];
except:
        gtest_mbdyn_args = ''
try:
        mbdyn_exec = os.environ['MBDYN_EXEC'];
except:
        mbdyn_exec = 'mbdyn'

os.system(mbdyn_exec + ' ' + gtest_mbdyn_args + ' -f cosim-output.mbd -o output -N 1 > output.txt 2>&1 &');

# wait for MBDyn to start and create sockets...
time.sleep(2);

# create input socket (4 double: Fx, Fy of 2 sensors)
# NOTE: what's out for MBDyn is in for this script
print('input socket: %s' % path_out);
s_in = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM);
s_in.connect(path_out);
s_in_bufsize = 6*8;


# step counter
i = 0
while 1:
        # read motion (x, x_prime) from body's output element
        buf_in = bytearray(s_in_bufsize)
        rc = s_in.recv_into(buf_in, s_in_bufsize)
        if (rc != s_in_bufsize):
                print("received %d bytes; closing..." % rc)
                break

        x1, x2, x3, theta1, theta2, theta3 = struct.unpack("dddddd", buf_in)

        print("i=%d x={%e,%e,%e} theta={%e,%e,%e} " % (i, x1, x2, x3, theta1, theta2, theta3))

        i = i + 1

s_in.close()
os.rmdir(tmpdir);
