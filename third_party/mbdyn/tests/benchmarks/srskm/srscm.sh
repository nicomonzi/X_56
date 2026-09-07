#!/bin/bash

FNAME=x
TMPFNAME=tmp
DT="4e-3"

awk -v DT="$DT" '$1 == 300 {
	printf "%+20.12e%+20.12e\n", DT*i++, $2;
}' "${FNAME}.mov" > "${TMPFNAME}.y.dat"

awk '$1 == 10 {
	phi = $16;
	while (phi - phi0 > 180) {
		phi -= 360;
	}
	while (phi - phi0 < -180) {
		phi += 360;
	}
	printf "%+20.12e\n", phi/57.2957795130823;
	phi0 = phi
}' "${FNAME}.jnt" > "${TMPFNAME}.theta.dat"

paste "${TMPFNAME}.y.dat" "${TMPFNAME}.theta.dat" > SpatialMBDyn.dat

rm -f "${TMPFNAME}.y.dat" "${TMPFNAME}.theta.dat"

