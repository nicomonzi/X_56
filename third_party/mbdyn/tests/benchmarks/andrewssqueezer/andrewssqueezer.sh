#!/bin/bash
#
# MBDyn (C) is a multibody analysis code. 
# http://www.mbdyn.org
# 
# Copyright (C) 1996-2021
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

if test "x$1" = x ; then
	echo "need <filename>.jnt"
	exit 1
fi
JNT=$1

if test `basename "$JNT" .jnt` == "$JNT" ; then
	echo "need <filename>.jnt"
	exit 1
fi

if ! test -f $JNT ; then
	echo "file $JNT does not exist"
	exit 1
fi

awk -v L_OF="0.007" -v every=250 \
'$1==0 {if ((i++%every) == 0) {printf "%20.12e%20.12e%20.12e\n", 1e-6*(i - 1), L_OF*cos($19), L_OF*sin($19)}}' \
$JNT > mov_F.txt
