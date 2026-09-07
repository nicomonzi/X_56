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

        body: 0, 0,
		M,
		reference, node, null,
		diag, J1, J1, J2;
	body: 1, 1,
		M,
		reference, node, null,
		diag, J1, J1, J2;
	body: 2, 2,
		M,
		reference, node, null,
		diag, J1, J1, J2;

	body: 101, 101,
		M,
		reference, node, null,
		diag, J2, J1, J1;
	body: 102, 102,
		M,
		reference, node, null,
		diag, J2, J1, J1;

	joint: 0, total pin joint,
		0,
			position, reference, global, 0., 0., 0.,
			position orientation, reference, global, eye,
			rotation orientation, reference, global, eye,
		position, reference, global, 0., 0., 0.,
		position orientation, reference, global, eye,
		rotation orientation, reference, global, eye,
		position constraint, 1, 1, 1, null,
		orientation constraint, 1, 0, 1, null;

	joint: 1, total pin joint,
		1,
			position, reference, global, L, 0., 0.,
			position orientation, reference, global, eye,
			rotation orientation, reference, global, eye,
		position, reference, global, L, 0., 0.,
		position orientation, reference, global, eye,
		rotation orientation, reference, global, eye,
		position constraint, 1, 0, 1, null,
		orientation constraint, 0, 0, 0, null;
	joint: 101, total joint,
		0,
			position, reference, global, 0., 0., L,
			position orientation, reference, global, eye,
			rotation orientation, reference, global, eye,
		101,
			position, reference, global, 0., 0., L,
			position orientation, reference, global, eye,
			rotation orientation, reference, global, eye,
		position constraint, 1, 1, 1, null,
		orientation constraint, 1, 0, 1, null;
	joint: 1101, total joint,
		1,
			position, reference, global, L, 0., L,
			position orientation, reference, global, eye,
			rotation orientation, reference, global, eye,
		101,
			position, reference, global, L, 0., L,
			position orientation, reference, global, eye,
			rotation orientation, reference, global, eye,
		position constraint, 1, 1, 1, null,
		orientation constraint, 1, 0, 1, null;

	joint: 2, total pin joint,
		2,
			position, reference, global, 2*L, 0., 0.,
			position orientation, reference, global, eye,
			rotation orientation, reference, global, eye,
		position, reference, global, 2*L, 0., 0.,
		position orientation, reference, global, eye,
		rotation orientation, reference, global, eye,
		position constraint, 1, 0, 1, null,
		orientation constraint, 0, 0, 0, null;
	joint: 102, total joint,
		1,
			position, reference, global, L, 0., L,
			position orientation, reference, global, eye,
			rotation orientation, reference, global, eye,
		102,
			position, reference, global, L, 0., L,
			position orientation, reference, global, eye,
			rotation orientation, reference, global, eye,
		position constraint, 1, 1, 1, null,
		orientation constraint, 1, 0, 1, null;
	joint: 1102, total joint,
		2,
			position, reference, global, 2*L, 0., L,
			position orientation, reference, global, eye,
			rotation orientation, reference, global, eye,
		102,
			position, reference, global, 2*L, 0., L,
			position orientation, reference, global, eye,
			rotation orientation, reference, global, eye,
		position constraint, 1, 1, 1, null,
		orientation constraint, 1, 0, 1, null;

	gravity: 0., 0., -1., 9.81;