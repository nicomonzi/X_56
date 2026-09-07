# MBDyn (C) is a multibody analysis code. 
# http://www.mbdyn.org
# 
# Copyright (C) 1996-2025
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

        joint: 0, total pin joint,
		0,
			position, reference, 1, null,
			position orientation, reference, 1, eye,
			rotation orientation, reference, 1, eye,
		position, reference, 1, null,
		position orientation, reference, 1, eye,
		rotation orientation, reference, 1, eye,
		position constraint, 1, 1, 1, null,
		orientation constraint, 1, 1, 1, null;

	include: "princeton.elm";

	force: N_BEAM_NODES, absolute,
		N_BEAM_NODES,
			position, reference, node, null,
		#reference, global,
			# 0., 0., -1.,
			0., sin(THETA*deg2rad), cos(THETA*deg2rad),
			cosine, 0., pi/.5, P/2, half, 0.;