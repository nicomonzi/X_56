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

        joint: 999, clamp, 999, node, node;

	driven: 998,
		string, "Time < TG/4",
		# string, "(Time > TG+.01) && (Time < TG/4)",
	joint: 998, rod,
		999,
			position, reference, SHAFT_R, L/2, 0., -1.,
		N_BEAMS,
			position, reference, SHAFT_R, L/2, 0., 0.,
		from nodes,
		linear viscous isotropic, 1e4;
			

	joint: SHAFT_R, total joint,
		999,
			position, reference, SHAFT_R, null,
			position orientation, reference, SHAFT_R, eye,
			rotation orientation, reference, SHAFT_R, eye,
		SHAFT_R,
			position, reference, SHAFT_R, null,
			position orientation, reference, SHAFT_R, eye,
			rotation orientation, reference, SHAFT_R, eye,
		position constraint, 1, 1, 1, null,
		orientation constraint, angular velocity, 1, 1,
			1., 0., 0.,
				# const, 0.;
				array, 2,
					cosine, T0, pi/T1, A1*OMEGA/2, half, 0.,
					cosine, T2, pi/(T3 - T2), (A2 - A1)*OMEGA/2, half, 0.;

	joint: SHAFT_T, total joint,
		999,
			position, reference, SHAFT_T, null,
			position orientation, reference, SHAFT_T, eye,
			rotation orientation, reference, SHAFT_T, eye,
		SHAFT_T,
			position, reference, SHAFT_T, null,
			position orientation, reference, SHAFT_T, eye,
			rotation orientation, reference, SHAFT_T, eye,
		position constraint, 0, 1, 1, null,
		orientation constraint, 0, 1, 1, null;

	include: "rotatingshaft_m.elm";
	
	gravity: 0., 0., -1,
		# const, 9.81;
		cosine, TG, pi/(-TG/2), 9.81/2, half, 0.;