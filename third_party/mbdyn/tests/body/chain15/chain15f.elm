# $Header: /var/cvs/mbdyn/mbdyn/mbdyn-1.0/tests/body/chain15/chain15f.elm,v 1.7 2017/01/12 15:00:35 masarati Exp $
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

# Author: Pierangelo Masarati <masarati@aero.polimi.it>
# From: C. G. Franchi
# "A Highly Redundant Coordinate Formulation for Constrained Rigid Bodies"
# Meccanica 30: 17-35, 1995.                                     1

joint: 1005, distance,
	0,
		position, reference, global, 0., DY, DELTA,
	5,
		position, reference, node, 0., 0., -L/2.,
	cosine, .01*T, pi/(.8*T), -sqrt(DY^2 + (5.*L)^2)/2, half,
		sqrt(DY^2 + (5.*L + DELTA)^2);

joint: 1010, distance,
	0,
		position, reference, global, 0., 2*DY, DELTA,
	10,
		position, reference, node, 0., 0., -L/2.,
	cosine, .01*T, pi/(.8*T), -sqrt((2*DY)^2 + (10.*L)^2)/2, half,
		sqrt((2*DY)^2 + (10.*L + DELTA)^2);

joint: 1015, distance,
	0,
		position, reference, global, 0., 3*DY, DELTA,
	15,
		position, reference, node, 0., 0., -L/2.,
	cosine, .01*T, pi/(.8*T), -sqrt((3*DY)^2 + (15.*L)^2)/2, half,
		sqrt((3*DY)^2 + (15.*L + DELTA)^2);

# vim:ft=mbd
