#!/bin/sh
#
# $Header: /var/cvs/mbdyn/mbdyn/mbdyn-1.0/tests/plates/sze2004/sze35/sze35_d.sh,v 1.1 2010/09/13 23:40:29 masarati Exp $

if test "x$1" = "x" ; then
	echo "usage: $0 <file without extension>"
	exit 1
fi

MOV="$1.mov"
FRC="$1.jnt"
OUT="`basename $1`.dat"

if test ! -f $MOV ; then
	echo "need $MOV file"
	exit 1
fi

if test ! -f $FRC ; then
	echo "need $FRC file"
	exit 1
fi

# POINT_A = 999998;
# POINT_B = 999999;
# POINT_C = 0;

awk '$1==0 {uC=4.953-$2} $1==999998 {wA=$4-4.953} $1==999999 {uB=4.953-$2; printf "%.6e %.6e %.6e\n", wA, uB, uC}' $MOV > tmp.mov
awk '$1==999998 {printf "%.6e\n", -$4}' $FRC > tmp.frc
paste tmp.frc tmp.mov > $OUT
rm -f tmp.frc tmp.mov

