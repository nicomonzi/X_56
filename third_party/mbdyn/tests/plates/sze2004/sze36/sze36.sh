#!/bin/sh
#
# $Header: /var/cvs/mbdyn/mbdyn/mbdyn-1.0/tests/plates/sze2004/sze36/sze36.sh,v 1.1 2010/09/12 07:37:33 masarati Exp $

if test "x$1" = "x" ; then
	echo "usage: $0 <file without extension>"
	exit 1
fi

MOV="$1.mov"
FRC="$1.frc"
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

awk '$1==999998 {wA=100-$4} $1==999999 {uB=$2-100; printf "%.6e %.6e\n", wA, uB}' $MOV > tmp.mov
awk '$1==999998 {printf "%.6e\n", -$5}' $FRC > tmp.frc
paste tmp.frc tmp.mov > $OUT
rm -f tmp.frc tmp.mov

