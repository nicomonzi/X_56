#!/bin/sh

# $Header: /var/cvs/mbdyn/mbdyn/mbdyn-1.0/tests/invdyn/pa10_traj.sh,v 1.6 2017/01/12 15:03:27 masarati Exp $

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

if test "x$1" = "x" ; then
	echo "usage: $0 [--delta DELTA] RESULTS_FILE"
	exit 1
fi

DELTA=1

while test "x$1" != "x" ; do
	case $1 in
	--delta)
		shift
		if test "x$1" = "x" ; then
			echo "usage: $0 [--delta DELTA] RESULTS_FILE"
			exit 1
		fi
		if test ! $1 -gt 0 ; then
			echo "error: DELTA must be greater than 0"
			exit 1
		fi
		DELTA=$1
		;;

	*)
		IN="$1.mov"
		if test ! -f "$IN" ; then
			echo "error: $IN must be a valid file"
			exit 1
		fi
		break
		;;
	esac
	shift
done

BASE=pa10.awk
OUT=pa10_traj.awk

cp -f $BASE $OUT

awk -v DELTA=$DELTA ' \
	BEGIN { \
		print "BEGIN {"; \
		print "\tedgeprop_add(\"traj\", 14, 1);"; \
		print ""; \
	} \
	$1 == 700 && (i++)%DELTA == 0 { \
		printf "\tnode_add(\"traj_%d\", 0, %e, %e, %e, \"hide\");\n", i, $2, $3, $4 - 0.005; \
	} \
	END { \
		print ""; \
		for (j = 1; j < i; j += DELTA) { \
			printf "\tedge_add(\"traj_%d\", \"traj_%d\", \"traj_%d\", \"traj\");\n", j, j, j+DELTA; \
		} \
		print "}"; \
	}' $IN >> $OUT
