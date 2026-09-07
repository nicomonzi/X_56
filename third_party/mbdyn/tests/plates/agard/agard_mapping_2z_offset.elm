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
# AGARD 445.6 wing model
# Author: Pierangelo Masarati
#         based on the Code Aster model provided by Giulio Romanelli

force: 111, external structural mapping,
	socket,
	create, yes,
	path, "/tmp/mbdyn_agard.sock",
	no signal,
	coupling,
		# loose,
		tight,
	reference node, 1,
	use reference node forces, no,
	points number, 242,
		1, offset, 0., 0., -.01, offset, 0., 0., .01,
		2, offset, 0., 0., -.01, offset, 0., 0., .01,
		3, offset, 0., 0., -.01, offset, 0., 0., .01,
		4, offset, 0., 0., -.01, offset, 0., 0., .01,
		5, offset, 0., 0., -.01, offset, 0., 0., .01,
		6, offset, 0., 0., -.01, offset, 0., 0., .01,
		7, offset, 0., 0., -.01, offset, 0., 0., .01,
		8, offset, 0., 0., -.01, offset, 0., 0., .01,
		9, offset, 0., 0., -.01, offset, 0., 0., .01,
		10, offset, 0., 0., -.01, offset, 0., 0., .01,
		11, offset, 0., 0., -.01, offset, 0., 0., .01,
		12, offset, 0., 0., -.01, offset, 0., 0., .01,
		13, offset, 0., 0., -.01, offset, 0., 0., .01,
		14, offset, 0., 0., -.01, offset, 0., 0., .01,
		15, offset, 0., 0., -.01, offset, 0., 0., .01,
		16, offset, 0., 0., -.01, offset, 0., 0., .01,
		17, offset, 0., 0., -.01, offset, 0., 0., .01,
		18, offset, 0., 0., -.01, offset, 0., 0., .01,
		19, offset, 0., 0., -.01, offset, 0., 0., .01,
		20, offset, 0., 0., -.01, offset, 0., 0., .01,
		21, offset, 0., 0., -.01, offset, 0., 0., .01,
		22, offset, 0., 0., -.01, offset, 0., 0., .01,
		23, offset, 0., 0., -.01, offset, 0., 0., .01,
		24, offset, 0., 0., -.01, offset, 0., 0., .01,
		25, offset, 0., 0., -.01, offset, 0., 0., .01,
		26, offset, 0., 0., -.01, offset, 0., 0., .01,
		27, offset, 0., 0., -.01, offset, 0., 0., .01,
		28, offset, 0., 0., -.01, offset, 0., 0., .01,
		29, offset, 0., 0., -.01, offset, 0., 0., .01,
		30, offset, 0., 0., -.01, offset, 0., 0., .01,
		31, offset, 0., 0., -.01, offset, 0., 0., .01,
		32, offset, 0., 0., -.01, offset, 0., 0., .01,
		33, offset, 0., 0., -.01, offset, 0., 0., .01,
		34, offset, 0., 0., -.01, offset, 0., 0., .01,
		35, offset, 0., 0., -.01, offset, 0., 0., .01,
		36, offset, 0., 0., -.01, offset, 0., 0., .01,
		37, offset, 0., 0., -.01, offset, 0., 0., .01,
		38, offset, 0., 0., -.01, offset, 0., 0., .01,
		39, offset, 0., 0., -.01, offset, 0., 0., .01,
		40, offset, 0., 0., -.01, offset, 0., 0., .01,
		41, offset, 0., 0., -.01, offset, 0., 0., .01,
		42, offset, 0., 0., -.01, offset, 0., 0., .01,
		43, offset, 0., 0., -.01, offset, 0., 0., .01,
		44, offset, 0., 0., -.01, offset, 0., 0., .01,
		45, offset, 0., 0., -.01, offset, 0., 0., .01,
		46, offset, 0., 0., -.01, offset, 0., 0., .01,
		47, offset, 0., 0., -.01, offset, 0., 0., .01,
		48, offset, 0., 0., -.01, offset, 0., 0., .01,
		49, offset, 0., 0., -.01, offset, 0., 0., .01,
		50, offset, 0., 0., -.01, offset, 0., 0., .01,
		51, offset, 0., 0., -.01, offset, 0., 0., .01,
		52, offset, 0., 0., -.01, offset, 0., 0., .01,
		53, offset, 0., 0., -.01, offset, 0., 0., .01,
		54, offset, 0., 0., -.01, offset, 0., 0., .01,
		55, offset, 0., 0., -.01, offset, 0., 0., .01,
		56, offset, 0., 0., -.01, offset, 0., 0., .01,
		57, offset, 0., 0., -.01, offset, 0., 0., .01,
		58, offset, 0., 0., -.01, offset, 0., 0., .01,
		59, offset, 0., 0., -.01, offset, 0., 0., .01,
		60, offset, 0., 0., -.01, offset, 0., 0., .01,
		61, offset, 0., 0., -.01, offset, 0., 0., .01,
		62, offset, 0., 0., -.01, offset, 0., 0., .01,
		63, offset, 0., 0., -.01, offset, 0., 0., .01,
		64, offset, 0., 0., -.01, offset, 0., 0., .01,
		65, offset, 0., 0., -.01, offset, 0., 0., .01,
		66, offset, 0., 0., -.01, offset, 0., 0., .01,
		67, offset, 0., 0., -.01, offset, 0., 0., .01,
		68, offset, 0., 0., -.01, offset, 0., 0., .01,
		69, offset, 0., 0., -.01, offset, 0., 0., .01,
		70, offset, 0., 0., -.01, offset, 0., 0., .01,
		71, offset, 0., 0., -.01, offset, 0., 0., .01,
		72, offset, 0., 0., -.01, offset, 0., 0., .01,
		73, offset, 0., 0., -.01, offset, 0., 0., .01,
		74, offset, 0., 0., -.01, offset, 0., 0., .01,
		75, offset, 0., 0., -.01, offset, 0., 0., .01,
		76, offset, 0., 0., -.01, offset, 0., 0., .01,
		77, offset, 0., 0., -.01, offset, 0., 0., .01,
		78, offset, 0., 0., -.01, offset, 0., 0., .01,
		79, offset, 0., 0., -.01, offset, 0., 0., .01,
		80, offset, 0., 0., -.01, offset, 0., 0., .01,
		81, offset, 0., 0., -.01, offset, 0., 0., .01,
		82, offset, 0., 0., -.01, offset, 0., 0., .01,
		83, offset, 0., 0., -.01, offset, 0., 0., .01,
		84, offset, 0., 0., -.01, offset, 0., 0., .01,
		85, offset, 0., 0., -.01, offset, 0., 0., .01,
		86, offset, 0., 0., -.01, offset, 0., 0., .01,
		87, offset, 0., 0., -.01, offset, 0., 0., .01,
		88, offset, 0., 0., -.01, offset, 0., 0., .01,
		89, offset, 0., 0., -.01, offset, 0., 0., .01,
		90, offset, 0., 0., -.01, offset, 0., 0., .01,
		91, offset, 0., 0., -.01, offset, 0., 0., .01,
		92, offset, 0., 0., -.01, offset, 0., 0., .01,
		93, offset, 0., 0., -.01, offset, 0., 0., .01,
		94, offset, 0., 0., -.01, offset, 0., 0., .01,
		95, offset, 0., 0., -.01, offset, 0., 0., .01,
		96, offset, 0., 0., -.01, offset, 0., 0., .01,
		97, offset, 0., 0., -.01, offset, 0., 0., .01,
		98, offset, 0., 0., -.01, offset, 0., 0., .01,
		99, offset, 0., 0., -.01, offset, 0., 0., .01,
		100, offset, 0., 0., -.01, offset, 0., 0., .01,
		101, offset, 0., 0., -.01, offset, 0., 0., .01,
		102, offset, 0., 0., -.01, offset, 0., 0., .01,
		103, offset, 0., 0., -.01, offset, 0., 0., .01,
		104, offset, 0., 0., -.01, offset, 0., 0., .01,
		105, offset, 0., 0., -.01, offset, 0., 0., .01,
		106, offset, 0., 0., -.01, offset, 0., 0., .01,
		107, offset, 0., 0., -.01, offset, 0., 0., .01,
		108, offset, 0., 0., -.01, offset, 0., 0., .01,
		109, offset, 0., 0., -.01, offset, 0., 0., .01,
		110, offset, 0., 0., -.01, offset, 0., 0., .01,
		111, offset, 0., 0., -.01, offset, 0., 0., .01,
		112, offset, 0., 0., -.01, offset, 0., 0., .01,
		113, offset, 0., 0., -.01, offset, 0., 0., .01,
		114, offset, 0., 0., -.01, offset, 0., 0., .01,
		115, offset, 0., 0., -.01, offset, 0., 0., .01,
		116, offset, 0., 0., -.01, offset, 0., 0., .01,
		117, offset, 0., 0., -.01, offset, 0., 0., .01,
		118, offset, 0., 0., -.01, offset, 0., 0., .01,
		119, offset, 0., 0., -.01, offset, 0., 0., .01,
		120, offset, 0., 0., -.01, offset, 0., 0., .01,
		121, offset, 0., 0., -.01, offset, 0., 0., .01,
	# echo, "agard_mb_2z.dat", precision, 16, surface, "agard_solid.dat", output, "agard_H_2z_offset.dat", order, 2, basenode, 12, weight, 2, stop;
	# echo, "agard_mb_2z.dat", precision, 16, surface, "agardSurfaceNew.dat", output, "agard_H_2z_offset.dat", order, 2, basenode, 12, weight, 2, stop;
	mapped points number, from file,
	sparse mapping file, "agard_H_2z_offset.dat";

# vim:ft=mbd

