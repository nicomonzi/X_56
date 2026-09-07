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

import socket
import struct

# create input socket (2 double: x, x_prime)
# s_in = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM);
# s_in.connect("./mbdyn.body.sock");
s_in = socket.socket(socket.AF_INET, socket.SOCK_DGRAM);
IN_PORT = 10011 
s_in.bind(("127.0.0.1", IN_PORT))
print(f"udpstream.py: UDP server up and listening on port {IN_PORT}")
s_in_bufsize = 2*8

# create output socket (1 double: f)
# s_out = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM);
# s_out.connect("./mbdyn.spring.sock");
OUT_PORT = 8005
s_out = socket.socket(socket.AF_INET, socket.SOCK_DGRAM);
s_out_bufsize = 8

# gains
KP = 1.e2
KD = 2.e-1

# initial value
f = 0.

# step counter
i = 0

# Pedantic output
DEBUG = True

while 1:
    # send force to body's file driver
        buf_out = bytearray(struct.pack("d", f))
        # rc = s_out.send(buf_out, s_out_bufsize)
        if DEBUG:
            print(f"sending {s_out_bufsize} bytes to localhost on port {OUT_PORT}")
        rc = s_out.sendto(buf_out, s_out_bufsize, ("127.0.0.1", OUT_PORT))	
        if (rc != s_out_bufsize):
            print(f"expecting to send {s_out_bufsize} bytes; closing...")
            break

        # read motion (x, x_prime) from body's output element
        buf_in = bytearray(s_in_bufsize)
        rc = s_in.recv_into(buf_in, s_in_bufsize)
        if DEBUG:
            print(f"received {rc} bytes to localhost on port {IN_PORT}")
        if (rc != s_in_bufsize):
            print(f"expecting {s_in_bufize} bytes; closing...")
            break

        x, xp = struct.unpack("dd", buf_in)

        # compute force for next step
        f = -KP*x - KD*xp

        # print("i=%d x=%e xp=%e f=%e" % (i, x, xp, f))

        i = i + 1

s_out.close()
s_in.close()

