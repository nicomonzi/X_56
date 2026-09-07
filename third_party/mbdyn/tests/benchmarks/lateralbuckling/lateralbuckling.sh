#!/bin/bash

# NOTE: must be odd!
BEAM_NBEAMS=5
LINK_NBEAMS=1
CRANK_NBEAMS=1
DT=0.0005

INFILE=zzzz
OUTFILE=MBDyn

# BeamMidAngVel
awk -v BEAM_NBEAMS=$BEAM_NBEAMS -v DT=$DT '
$1==1000 + BEAM_NBEAMS {
	printf "%20.15f %20.15f\n", DT*i++, $17;
}' $INFILE.mov > ${OUTFILE}_BeamMidAngVel.dat


# BeamMidDisplacement
awk -v BEAM_NBEAMS=$BEAM_NBEAMS -v DT=$DT '
$1==1000 + BEAM_NBEAMS {
	printf "%20.15f %20.15f\n", DT*i++, $3;
}' $INFILE.mov > ${OUTFILE}_BeamMidDisplacement.dat

# BeamMidForce
awk -v BEAM_NBEAMS=$BEAM_NBEAMS -v DT=$DT '
$1==1000 + int(BEAM_NBEAMS/2) + 1 {
	printf "%20.15f %20.15f\n", DT*i++, ($4 + $10)/2;
}' $INFILE.act > ${OUTFILE}_BeamMidForce.dat


# BeamMidMoment
awk -v BEAM_NBEAMS=$BEAM_NBEAMS -v DT=$DT '
$1==1000 + int(BEAM_NBEAMS/2) + 1 {
	printf "%20.15f %20.15f\n", DT*i++, ($6 + $12)/2;
}' $INFILE.act > ${OUTFILE}_BeamMidMoment.dat


# BeamMidRotation
# from http://www.dymoresolutions.com/SectionBuilder/UsersManual/GeometricObjects/Orientations.pdf
# phi3 = atan2(-R31, R33)
# phi2 = atan2(R32, R33*cos(phi3) - R31*sin(phi3))
awk -v BEAM_NBEAMS=$BEAM_NBEAMS -v DT=$DT '
$1==1000 + BEAM_NBEAMS {
	R31 = $11;
	R32 = $12;
	R33 = $13;
	phi3 = atan2(-R31, R33);
	phi2 = atan2(R32, R33*cos(phi3) - R31*sin(phi3));
	printf "%20.15f %20.15f\n", DT*i++, phi2*57.2957795130823;
}' $INFILE.mov > ${OUTFILE}_BeamMidRotation.dat


# LinkMidForce
awk -v DT=$DT '
$1==2000 + 1 {
	printf "%20.15f %20.15f\n", DT*i++, ($2 + $8)/2;
}' $INFILE.act > ${OUTFILE}_LinkMidForce.dat


# LinkMidMoment
awk -v DT=$DT '
$1==2000 + 1 {
	printf "%20.15f %20.15f\n", DT*i++, ($7 + $13)/2;
}' $INFILE.act > ${OUTFILE}_LinkMidMoment.dat

