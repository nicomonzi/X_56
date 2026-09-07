#! /bin/bash

FNAME="x"
DT="0.004"
STEP=5236

awk -v DT="$DT" -v STEP="$STEP" '$1==2010 {
	if (i == 0) {
		x0 = $2;
		y0 = $4;
		z0 = -$3;
	}
	if (++i > STEP) {
		print DT*(i - STEP - 1), $2 - x0, $4 - y0, -$3 - z0, 0., 0., 0.;
	}
}' "${FNAME}.mov" > SensorBeam2TipDisplacementsMBDyn.mdt

awk -v DT="$DT" -v STEP="$STEP" '$1==1003 {
	if (++i > STEP) {
		print DT*(i - STEP - 1), ($2 + $8)/2, ($3 + $9)/2, ($4 + $10)/2, ($5 + $11)/2, ($6 + $12)/2, ($7 + $13)/2;
	}
}' "${FNAME}.act" > SensorBeam1MidForcesMBDyn.mdt

awk -v DT="$DT" -v STEP="$STEP" '$1==1005 {
	if (++i > STEP) {
		print DT*(i - STEP - 1), $8, $10, -$9, $11, $13, -$12;
	}
}' "${FNAME}.mov" > SensorBeam1MidVelocitiesMBDyn.mdt

awk -v DT="$DT" -v STEP="$STEP" '$1==4 {
	if (++i > STEP) {
		print DT*(i - STEP - 1), 0., 0., 0., -$19, 0., 0.;
	}
}' "${FNAME}.jnt" > SensorRotationDMBDyn.mdt

