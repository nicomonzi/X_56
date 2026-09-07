#!/bin/bash

OCTAVE_EXEC=${OCTAVE_EXEC:-octave}
MBDYN_EXEC=${MBDYN_EXEC:-mbdyn}

${OCTAVE_EXEC} ${GTEST_OCTAVE_ARGS} -q princeton_gen.m $*

declare -a PP
PP[1]="4.448"
PP[2]="8.896"
PP[3]="13.345"

for P in 1 2 3 ; do
	DATFILENAME=`printf "princeton_p%1d.dat" $P`
	echo "# theta u_2 u_3 phi" > $DATFILENAME

	DATFILENAME2=`printf "PrincetonBeamMBDyn_P%1d.dat" $P`
	echo "" > $DATFILENAME2

	THETA=0
	DTHETA=5
	while test $THETA -le 90 ; do
		FILENAME=`printf "princeton_theta%02d_p%1d" $THETA $P`

		MBDYNVARS="real THETA=$THETA; real P=${PP[$P]}" \
			${MBDYN_EXEC} ${GTEST_MBDYN_ARGS} -f princeton -o $FILENAME -ss

		tail -1 $FILENAME.mov | awk '{print '$THETA', $3, $4, atan2(-$10, $13)}' >> $DATFILENAME

		tail -1 $FILENAME.mov | awk '{print '$THETA', $2, $3, $4, $5, $6, $7}' >> $DATFILENAME2

		THETA=`expr $THETA + $DTHETA`
	done
done


paste princeton_p*.dat | awk '{
	if ($1 == "#") {
		print "" > "PrincetonBeamMBDyn_u2.dat";
		print "" > "PrincetonBeamMBDyn_u3.dat";
		print "" > "PrincetonBeamMBDyn_phi.dat";
	} else {
		printf "%20.15f%20.15f%20.15f\n", $2, $6, $10 > "PrincetonBeamMBDyn_u2.dat"; 
		printf "%20.15f%20.15f%20.15f\n", $3, $7, $11 > "PrincetonBeamMBDyn_u3.dat"; 
		printf "%20.15f%20.15f%20.15f\n", $4, $8, $12 > "PrincetonBeamMBDyn_phi.dat"; 
	}
}'

gnuplot princeton.plt

rm -f princeton_theta*
