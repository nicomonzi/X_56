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
	points number, 121,
		1, offset, null,
		2, offset, null,
		3, offset, null,
		4, offset, null,
		5, offset, null,
		6, offset, null,
		7, offset, null,
		8, offset, null,
		9, offset, null,
		10, offset, null,
		11, offset, null,
		12, offset, null,
		13, offset, null,
		14, offset, null,
		15, offset, null,
		16, offset, null,
		17, offset, null,
		18, offset, null,
		19, offset, null,
		20, offset, null,
		21, offset, null,
		22, offset, null,
		23, offset, null,
		24, offset, null,
		25, offset, null,
		26, offset, null,
		27, offset, null,
		28, offset, null,
		29, offset, null,
		30, offset, null,
		31, offset, null,
		32, offset, null,
		33, offset, null,
		34, offset, null,
		35, offset, null,
		36, offset, null,
		37, offset, null,
		38, offset, null,
		39, offset, null,
		40, offset, null,
		41, offset, null,
		42, offset, null,
		43, offset, null,
		44, offset, null,
		45, offset, null,
		46, offset, null,
		47, offset, null,
		48, offset, null,
		49, offset, null,
		50, offset, null,
		51, offset, null,
		52, offset, null,
		53, offset, null,
		54, offset, null,
		55, offset, null,
		56, offset, null,
		57, offset, null,
		58, offset, null,
		59, offset, null,
		60, offset, null,
		61, offset, null,
		62, offset, null,
		63, offset, null,
		64, offset, null,
		65, offset, null,
		66, offset, null,
		67, offset, null,
		68, offset, null,
		69, offset, null,
		70, offset, null,
		71, offset, null,
		72, offset, null,
		73, offset, null,
		74, offset, null,
		75, offset, null,
		76, offset, null,
		77, offset, null,
		78, offset, null,
		79, offset, null,
		80, offset, null,
		81, offset, null,
		82, offset, null,
		83, offset, null,
		84, offset, null,
		85, offset, null,
		86, offset, null,
		87, offset, null,
		88, offset, null,
		89, offset, null,
		90, offset, null,
		91, offset, null,
		92, offset, null,
		93, offset, null,
		94, offset, null,
		95, offset, null,
		96, offset, null,
		97, offset, null,
		98, offset, null,
		99, offset, null,
		100, offset, null,
		101, offset, null,
		102, offset, null,
		103, offset, null,
		104, offset, null,
		105, offset, null,
		106, offset, null,
		107, offset, null,
		108, offset, null,
		109, offset, null,
		110, offset, null,
		111, offset, null,
		112, offset, null,
		113, offset, null,
		114, offset, null,
		115, offset, null,
		116, offset, null,
		117, offset, null,
		118, offset, null,
		119, offset, null,
		120, offset, null,
		121, offset, null,
	# echo, "agard_mb.dat", surface, "agard_solid.dat", output, "agard_H_no_offset.dat", order, 2, basenode, 12, weight, 2, stop;
	mapped points number, 840,
	sparse mapping file, "agard_H_no_offset.dat";

# vim:ft=mbd

