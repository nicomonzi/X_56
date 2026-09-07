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
		1, offset, null, offset, 0., 0., .02,
		2, offset, null, offset, 0., 0., .02,
		3, offset, null, offset, 0., 0., .02,
		4, offset, null, offset, 0., 0., .02,
		5, offset, null, offset, 0., 0., .02,
		6, offset, null, offset, 0., 0., .02,
		7, offset, null, offset, 0., 0., .02,
		8, offset, null, offset, 0., 0., .02,
		9, offset, null, offset, 0., 0., .02,
		10, offset, null, offset, 0., 0., .02,
		11, offset, null, offset, 0., 0., .02,
		12, offset, null, offset, 0., 0., .02,
		13, offset, null, offset, 0., 0., .02,
		14, offset, null, offset, 0., 0., .02,
		15, offset, null, offset, 0., 0., .02,
		16, offset, null, offset, 0., 0., .02,
		17, offset, null, offset, 0., 0., .02,
		18, offset, null, offset, 0., 0., .02,
		19, offset, null, offset, 0., 0., .02,
		20, offset, null, offset, 0., 0., .02,
		21, offset, null, offset, 0., 0., .02,
		22, offset, null, offset, 0., 0., .02,
		23, offset, null, offset, 0., 0., .02,
		24, offset, null, offset, 0., 0., .02,
		25, offset, null, offset, 0., 0., .02,
		26, offset, null, offset, 0., 0., .02,
		27, offset, null, offset, 0., 0., .02,
		28, offset, null, offset, 0., 0., .02,
		29, offset, null, offset, 0., 0., .02,
		30, offset, null, offset, 0., 0., .02,
		31, offset, null, offset, 0., 0., .02,
		32, offset, null, offset, 0., 0., .02,
		33, offset, null, offset, 0., 0., .02,
		34, offset, null, offset, 0., 0., .02,
		35, offset, null, offset, 0., 0., .02,
		36, offset, null, offset, 0., 0., .02,
		37, offset, null, offset, 0., 0., .02,
		38, offset, null, offset, 0., 0., .02,
		39, offset, null, offset, 0., 0., .02,
		40, offset, null, offset, 0., 0., .02,
		41, offset, null, offset, 0., 0., .02,
		42, offset, null, offset, 0., 0., .02,
		43, offset, null, offset, 0., 0., .02,
		44, offset, null, offset, 0., 0., .02,
		45, offset, null, offset, 0., 0., .02,
		46, offset, null, offset, 0., 0., .02,
		47, offset, null, offset, 0., 0., .02,
		48, offset, null, offset, 0., 0., .02,
		49, offset, null, offset, 0., 0., .02,
		50, offset, null, offset, 0., 0., .02,
		51, offset, null, offset, 0., 0., .02,
		52, offset, null, offset, 0., 0., .02,
		53, offset, null, offset, 0., 0., .02,
		54, offset, null, offset, 0., 0., .02,
		55, offset, null, offset, 0., 0., .02,
		56, offset, null, offset, 0., 0., .02,
		57, offset, null, offset, 0., 0., .02,
		58, offset, null, offset, 0., 0., .02,
		59, offset, null, offset, 0., 0., .02,
		60, offset, null, offset, 0., 0., .02,
		61, offset, null, offset, 0., 0., .02,
		62, offset, null, offset, 0., 0., .02,
		63, offset, null, offset, 0., 0., .02,
		64, offset, null, offset, 0., 0., .02,
		65, offset, null, offset, 0., 0., .02,
		66, offset, null, offset, 0., 0., .02,
		67, offset, null, offset, 0., 0., .02,
		68, offset, null, offset, 0., 0., .02,
		69, offset, null, offset, 0., 0., .02,
		70, offset, null, offset, 0., 0., .02,
		71, offset, null, offset, 0., 0., .02,
		72, offset, null, offset, 0., 0., .02,
		73, offset, null, offset, 0., 0., .02,
		74, offset, null, offset, 0., 0., .02,
		75, offset, null, offset, 0., 0., .02,
		76, offset, null, offset, 0., 0., .02,
		77, offset, null, offset, 0., 0., .02,
		78, offset, null, offset, 0., 0., .02,
		79, offset, null, offset, 0., 0., .02,
		80, offset, null, offset, 0., 0., .02,
		81, offset, null, offset, 0., 0., .02,
		82, offset, null, offset, 0., 0., .02,
		83, offset, null, offset, 0., 0., .02,
		84, offset, null, offset, 0., 0., .02,
		85, offset, null, offset, 0., 0., .02,
		86, offset, null, offset, 0., 0., .02,
		87, offset, null, offset, 0., 0., .02,
		88, offset, null, offset, 0., 0., .02,
		89, offset, null, offset, 0., 0., .02,
		90, offset, null, offset, 0., 0., .02,
		91, offset, null, offset, 0., 0., .02,
		92, offset, null, offset, 0., 0., .02,
		93, offset, null, offset, 0., 0., .02,
		94, offset, null, offset, 0., 0., .02,
		95, offset, null, offset, 0., 0., .02,
		96, offset, null, offset, 0., 0., .02,
		97, offset, null, offset, 0., 0., .02,
		98, offset, null, offset, 0., 0., .02,
		99, offset, null, offset, 0., 0., .02,
		100, offset, null, offset, 0., 0., .02,
		101, offset, null, offset, 0., 0., .02,
		102, offset, null, offset, 0., 0., .02,
		103, offset, null, offset, 0., 0., .02,
		104, offset, null, offset, 0., 0., .02,
		105, offset, null, offset, 0., 0., .02,
		106, offset, null, offset, 0., 0., .02,
		107, offset, null, offset, 0., 0., .02,
		108, offset, null, offset, 0., 0., .02,
		109, offset, null, offset, 0., 0., .02,
		110, offset, null, offset, 0., 0., .02,
		111, offset, null, offset, 0., 0., .02,
		112, offset, null, offset, 0., 0., .02,
		113, offset, null, offset, 0., 0., .02,
		114, offset, null, offset, 0., 0., .02,
		115, offset, null, offset, 0., 0., .02,
		116, offset, null, offset, 0., 0., .02,
		117, offset, null, offset, 0., 0., .02,
		118, offset, null, offset, 0., 0., .02,
		119, offset, null, offset, 0., 0., .02,
		120, offset, null, offset, 0., 0., .02,
		121, offset, null, offset, 0., 0., .02,
	# echo, "agard_mb_1z.dat", surface, "agard_solid.dat", output, "agard_H_1z_offset.dat", order, 2, basenode, 12, weight, 2, stop;
	# echo, "agard_mb_1z.dat", surface, "agardSurfaceNew.dat", output, "agard_H_1z_offset.dat", order, 2, basenode, 12, weight, 2, stop;
	mapped points number, from file,
	sparse mapping file, "agard_H_1z_offset.dat";

# vim:ft=mbd

