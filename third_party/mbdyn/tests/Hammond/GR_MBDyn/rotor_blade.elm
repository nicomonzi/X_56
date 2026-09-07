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
#
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
# ---------------------------------------------------------------------------
#
# 4-blade rotor on flexible support, from 
# Hammond, C. E., 1974. "An application of Floquet theory to prediction of mechanical instability".
# Journal of the American Helicopter Society, 19(4), pp. 14–23. doi:10.4050/JAHS.19.14.
#
# Author: Gianni Cassoni
# Based on a generic rotor model by Pierangelo Masarati
#
# ------ GR MODEL: blade elements ------

body: CURR_BLADE, CURR_BLADE,
	M_BLADE,
	reference, CURR_BLADE + BLADE_CM, null,
	diag, J_BLADE_CM_X, J_BLADE_CM_Y, J_BLADE_CM_Z;

joint: CURR_BLADE + BLADE_FLAP_OFFSET, total joint,
	HUB,
		position, reference, CURR_BLADE + BLADE_FLAP_OFFSET, null,
		position orientation, reference, CURR_BLADE + BLADE_FLAP_OFFSET, eye,
		rotation orientation, reference, CURR_BLADE + BLADE_FLAP_OFFSET, eye,
	CURR_BLADE,
		position, reference, CURR_BLADE + BLADE_FLAP_OFFSET, null,
		position orientation, reference, CURR_BLADE + BLADE_FLAP_OFFSET, eye,
		rotation orientation, reference, CURR_BLADE + BLADE_FLAP_OFFSET, eye,
	position constraint, 1, 1, 1, null,
	orientation constraint, 1, 1, 0, null;
joint: CURR_BLADE + BLADE_FLAP_OFFSET + 1, deformable hinge,
	HUB,
		position, reference, CURR_BLADE + BLADE_FLAP_OFFSET, null,
		orientation, reference, CURR_BLADE + BLADE_FLAP_OFFSET, eye,
	CURR_BLADE,
		position, reference, CURR_BLADE + BLADE_FLAP_OFFSET, null,
		orientation, reference, CURR_BLADE + BLADE_FLAP_OFFSET, eye,
	reference, CURR_CL_BLADE_ROOT;

# vim:ft=mbd

