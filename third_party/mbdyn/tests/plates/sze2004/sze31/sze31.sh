#!/bin/sh
#
# $Header: /var/cvs/mbdyn/mbdyn/mbdyn-1.0/tests/plates/sze2004/sze31/sze31.sh,v 1.2 2010/09/11 03:09:14 masarati Exp $

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

awk '$1==199 {printf "%.6e %.6e\n", 10 - $2, $4}' $MOV > tmp.mov
awk '$1==100 {printf "%.6e\n", $5*2/4}' $FRC > tmp.frc
paste tmp.frc tmp.mov > $OUT
rm -f tmp.frc tmp.mov

