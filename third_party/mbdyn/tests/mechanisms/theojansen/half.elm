# $Header: /var/cvs/mbdyn/mbdyn/mbdyn-1.0/tests/mechanisms/theojansen/half.elm,v 1.9 2017/01/12 15:05:53 masarati Exp $
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
# Author: Pierangelo Masarati <masarati@aero.polimi.it>
# Model: Theo Jansen 's Mechanism 
# From:	<http://il.youtube.com/watch?v=-GgOn66knqA>
# Data from: <http://www.mechanicalspider.com/elements/jpegs/jansenlinkage.jpg>
#
# uses "main_assy_*.ref" computed from the solution of "main_assy"

body: CURR_BAY + CURR_SIDE + LINK_1, CURR_BAY + CURR_SIDE + LINK_1, 1.e-3, null, eye;
body: CURR_BAY + CURR_SIDE + LINK_2, CURR_BAY + CURR_SIDE + LINK_2, 1.e-3, null, eye;
body: CURR_BAY + CURR_SIDE + LINK_3, CURR_BAY + CURR_SIDE + LINK_3, 1.e-3, null, eye;
body: CURR_BAY + CURR_SIDE + LINK_4, CURR_BAY + CURR_SIDE + LINK_4, 1.e-3, null, eye;
body: CURR_BAY + CURR_SIDE + TOP_BODY, CURR_BAY + CURR_SIDE + TOP_BODY, 1.e-3, null, eye;
body: CURR_BAY + CURR_SIDE + BOTTOM_BODY, CURR_BAY + CURR_SIDE + BOTTOM_BODY, 1.e-3, null, eye;

joint: CURR_BAY + CURR_SIDE + BAR + 1, revolute hinge,
	BAR,
		position, reference, CURR_BAY + CURR_SIDE + POINT_ABCD, null,
		orientation, reference, CURR_BAY + CURR_SIDE + POINT_ABCD,
			3, 0., 1., 0.,
			1, 1., 0., 0.,
	CURR_BAY + CURR_SIDE + TOP_BODY,
		position, reference, CURR_BAY + CURR_SIDE + POINT_ABCD, null,
		orientation, reference, CURR_BAY + CURR_SIDE + POINT_ABCD,
			3, 0., 1., 0.,
			1, 1., 0., 0.;

joint: CURR_BAY + CURR_SIDE + LINK_1 + 1, revolute hinge,
	CURR_BAY + CURR_SIDE + TOP_BODY,
		position, reference, CURR_BAY + CURR_SIDE + POINT_DEF, null,
		orientation, reference, CURR_BAY + CURR_SIDE + POINT_DEF,
			3, 0., 1., 0.,
			1, 1., 0., 0.,
	CURR_BAY + CURR_SIDE + LINK_1,
		position, reference, CURR_BAY + CURR_SIDE + POINT_DEF, null,
		orientation, reference, CURR_BAY + CURR_SIDE + POINT_DEF,
			3, 0., 1., 0.,
			1, 1., 0., 0.;
joint: CURR_BAY + CURR_SIDE + LINK_1 + 2, revolute hinge,
	CURR_BAY + CURR_SIDE + BOTTOM_BODY,
		position, reference, CURR_BAY + CURR_SIDE + POINT_FGH, null,
		orientation, reference, CURR_BAY + CURR_SIDE + POINT_FGH,
			3, 0., 1., 0.,
			1, 1., 0., 0.,
	CURR_BAY + CURR_SIDE + LINK_1,
		position, reference, CURR_BAY + CURR_SIDE + POINT_FGH, null,
		orientation, reference, CURR_BAY + CURR_SIDE + POINT_FGH,
			3, 0., 1., 0.,
			1, 1., 0., 0.;

joint: CURR_BAY + CURR_SIDE + LINK_2 + 1, spherical hinge,
	CURR_BAY + CURR_SIDE + TOP_BODY,
		position, reference, CURR_BAY + CURR_SIDE + POINT_ABCD, null,
		orientation, reference, CURR_BAY + CURR_SIDE + POINT_ABCD, eye,
	CURR_BAY + CURR_SIDE + LINK_2,
		position, reference, CURR_BAY + CURR_SIDE + POINT_ABCD, null,
		orientation, reference, CURR_BAY + CURR_SIDE + POINT_ABCD, eye;
joint: CURR_BAY + CURR_SIDE + LINK_2 + 2, spherical hinge,
	CURR_BAY + CURR_SIDE + BOTTOM_BODY,
		position, reference, CURR_BAY + CURR_SIDE + POINT_CGIK, null,
		orientation, reference, CURR_BAY + CURR_SIDE + POINT_CGIK, eye,
	CURR_BAY + CURR_SIDE + LINK_2,
		position, reference, CURR_BAY + CURR_SIDE + POINT_CGIK, null,
		orientation, reference, CURR_BAY + CURR_SIDE + POINT_CGIK, eye;

joint: CURR_BAY + CURR_SIDE + LINK_3 + 1, spherical hinge,
	CRANK,
		position, reference, CURR_BAY + POINT_JKM, null,
		orientation, reference, CURR_BAY + POINT_JKM, eye,
	CURR_BAY + CURR_SIDE + LINK_3,
		position, reference, CURR_BAY + POINT_JKM, null,
		orientation, reference, CURR_BAY + POINT_JKM, eye;

joint: CURR_BAY + CURR_SIDE + LINK_3 + 2, spherical hinge,
	CURR_BAY + CURR_SIDE + TOP_BODY,
		position, reference, CURR_BAY + CURR_SIDE + POINT_BEJ, null,
		orientation, reference, CURR_BAY + CURR_SIDE + POINT_BEJ, eye,
	CURR_BAY + CURR_SIDE + LINK_3,
		position, reference, CURR_BAY + CURR_SIDE + POINT_BEJ, null,
		orientation, reference, CURR_BAY + CURR_SIDE + POINT_BEJ, eye;

joint: CURR_BAY + CURR_SIDE + LINK_4 + 1, spherical hinge,
	CRANK,
		position, reference, CURR_BAY + POINT_JKM, null,
		orientation, reference, CURR_BAY + POINT_JKM, eye,
	CURR_BAY + CURR_SIDE + LINK_4,
		position, reference, CURR_BAY + POINT_JKM, null,
		orientation, reference, CURR_BAY + POINT_JKM, eye;
joint: CURR_BAY + CURR_SIDE + LINK_4 + 2, spherical hinge,
	CURR_BAY + CURR_SIDE + BOTTOM_BODY,
		position, reference, CURR_BAY + CURR_SIDE + POINT_CGIK, null,
		orientation, reference, CURR_BAY + CURR_SIDE + POINT_CGIK, eye,
	CURR_BAY + CURR_SIDE + LINK_4,
		position, reference, CURR_BAY + CURR_SIDE + POINT_CGIK, null,
		orientation, reference, CURR_BAY + CURR_SIDE + POINT_CGIK, eye;

# vim:ft=mbd
