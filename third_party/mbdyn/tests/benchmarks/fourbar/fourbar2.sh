#! /bin/bash

FNAME="x"
DT="0.004"
ISTEP=5236
DSTEP=2618
FSTEP=`expr $ISTEP + $DSTEP`
NBEAMS=5
BAR_1_MID_NODE=`expr 1000 + $NBEAMS`
BAR_1_MID_BEAM=`expr 1000 + $NBEAMS / 2 + 1`
POINT_C=3000
BAR_2_END_NODE=`expr 2000 + 2 \* $NBEAMS`
echo "BAR_2_END_NODE=$BAR_2_END_NODE"

awk -v DT="$DT" -v ISTEP="$ISTEP" -v FSTEP="$FSTEP" '$1==1 {
	if (i == 0) {
		print "# Time[s]";
	}
	if (++i > ISTEP && i <= FSTEP) {
		printf "%20.15f\n", DT*(i - ISTEP - 1);
	}
}' "${FNAME}.jnt" > Time_MBDyn.dat

awk -v ISTEP="$ISTEP" -v FSTEP="$FSTEP" '$1==1 {
	if (i == 0) {
		print "# CrankAngle[deg]";
	}
	if (++i > ISTEP && i <= FSTEP) {
		phi = -$19*57.2957795130823;
		while (phi - phi0 < -180) {
			phi += 360;
		}
		while (phi - phi0 > 180) {
			phi -= 360;
		}
		printf "%20.15f\n", phi;
		phi0 = phi;
	}
}' "${FNAME}.jnt" > CrankAngle_MBDyn.dat

awk -v ISTEP="$ISTEP" -v FSTEP="$FSTEP" -v BAR_1_MID_NODE="$BAR_1_MID_NODE" '$1==BAR_1_MID_NODE {
	if (i == 0) {
		print "# AngularVelocity[rad/s]";
	}
	if (++i > ISTEP && i <= FSTEP) {
		printf "%20.15f\n", -$12;
	}
}' "${FNAME}.mov" > Bar1AngularVelocity_MBDyn.dat

awk -v ISTEP="$ISTEP" -v FSTEP="$FSTEP" -v BAR_1_MID_NODE="$BAR_1_MID_NODE" '$1==BAR_1_MID_NODE {
	if (i == 0) {
		print "# Velocity[m/s]";
	}
	if (++i > ISTEP && i <= FSTEP) {
		printf "%20.15f\n", $8;
	}
}' "${FNAME}.mov" > Bar1MidVelocity_MBDyn.dat

awk -v ISTEP="$ISTEP" -v FSTEP="$FSTEP" -v BAR_1_MID_BEAM="$BAR_1_MID_BEAM" '$1==BAR_1_MID_BEAM {
	if (i == 0) {
		print "# Force[N]";
	}
	if (++i > ISTEP && i <= FSTEP) {
		printf "%20.15f\n", ($2 + $8)/2;
	}
}' "${FNAME}.act" > Bar1MidForce_MBDyn.dat

awk -v ISTEP="$ISTEP" -v FSTEP="$FSTEP" -v BAR_1_MID_BEAM="$BAR_1_MID_BEAM" '$1==BAR_1_MID_BEAM {
	if (i == 0) {
		print "# Moment[Nm]";
	}
	if (++i > ISTEP && i <= FSTEP) {
		printf "%20.15f\n", ($6 + $12)/2;
	}
}' "${FNAME}.act" > Bar1MidMoment_MBDyn.dat

awk -v ISTEP="$ISTEP" -v FSTEP="$FSTEP" -v POINT_C="$POINT_C" '$1==POINT_C {
	if (i == 0) {
		x0 = $2;
		print "# Displacement[m]";
	}
	if (++i > ISTEP && i <= FSTEP) {
		printf "%20.15f\n", $2 - x0;
	}
}' "${FNAME}.mov" > DispC_MBDyn.dat

awk -v ISTEP="$ISTEP" -v FSTEP="$FSTEP" '$1==4 {
	if (i == 0) {
		print "# Rotation[rad]";
	}
	if (++i > ISTEP && i <= FSTEP) {
		phi = -$19*57.2957795130823;
		while (phi - phi0 < -180) {
			phi += 360;
		}
		while (phi - phi0 > 180) {
			phi -= 360;
		}
		printf "%20.15f\n", phi/57.2957795130823;
		phi0 = phi;
	}
}' "${FNAME}.jnt" > RotationD_MBDyn.dat

#awk -v ISTEP="$ISTEP" -v FSTEP="$FSTEP" -v BAR_2_END_NODE="$BAR_2_END_NODE" '$1==BAR_2_END_NODE {
#	if (i == 0) {
#		print "# Euler-3-1-2-psi[deg] Euler-3-1-2-theta[deg] Euler-3-1-2-phi[deg]";
#	}
#	if (++i > ISTEP && i <= FSTEP) {
#		phi1 = $5;
#		phi2 = $6;
#		phi3 = $7;
#		dPhi = sqrt(phi1*phi1 + phi2*phi2 + phi3*phi3);
#		da = sin(dPhi)/dPhi;
#		db = (1 - cos(dPhi))/(dPhi*dPhi);
#		# psi = atan2(-R31, R33); Bauchau, 4.11.4 Euler angles: sequence 3-1-2 (4.78) pag. 140
#
#		R11 = 1. - db*(phi2*phi2 + phi3*phi3);
#		R21 = da*phi3 + db*phi2*phi1;
#		R31 = -da*phi2 + db*phi3*phi1;
#		R12 = -da*phi3 + db*phi1*phi2;
#		R22 = 1. - db*(phi3*phi3 + phi1*phi1);
#		R32 = da*phi1 + db*phi3*phi2;
#		R13 = da*phi2 + db*phi1*phi3;
#		R23 = -da*phi1 + db*phi2*phi3;
#		R33 = 1. - db*(phi1*phi1 + phi2*phi2);
#
#		psi = atan2(-R31, R33);
#		theta = atan2(R32, -R31*sin(psi) + R33*cos(psi));
#		phi = atan2(-R12, R22);
#		
#		printf "%20.15f%20.15f%20.15f\n", psi*57.2957795130823, theta*57.2957795130823, phi*57.2957795130823;
#	}
#}' "${FNAME}.mov" > Bar2RotC_MBDyn.dat

awk -v ISTEP="$ISTEP" -v FSTEP="$FSTEP" -v BAR_2_END_NODE="$BAR_2_END_NODE" '$1==BAR_2_END_NODE {
	if (i == 0) {
		print "# Euler-3-1-2-psi[deg] Euler-3-1-2-theta[deg] Euler-3-1-2-phi[deg]";
	}
	if (++i > ISTEP && i <= FSTEP) {
		R11 = $5;
		R21 = $8;
		R31 = $11;
		R12 = $6;
		R22 = $9;
		R32 = $12;
		R13 = $7;
		R23 = $10;
		R33 = $13;

		psi = atan2(-R31, R33);
		theta = atan2(R32, -R31*sin(psi) + R33*cos(psi));
		phi = atan2(-R12, R22);
		
		printf "%20.15f\n", phi*57.2957795130823;
	}
}' "${FNAME}.mov" > Bar2RotC_MBDyn.dat

awk -v ISTEP="$ISTEP" -v FSTEP="$FSTEP" -v BAR_2_END_NODE="$BAR_2_END_NODE" '$1==BAR_2_END_NODE {
	if (i == 0) {
		print "# Euler-3-1-2-psi[deg] Euler-3-1-2-theta[deg] Euler-3-1-2-phi[deg]";
	}
	if (++i > ISTEP && i <= FSTEP) {
		R11 = $5;
		R21 = $8;
		R31 = $11;
		R12 = $6;
		R22 = $9;
		R32 = $12;
		R13 = $7;
		R23 = $10;
		R33 = $13;

		psi = atan2(-R31, R33);
		theta = atan2(R32, -R31*sin(psi) + R33*cos(psi));
		phi = atan2(-R12, R22);
		
		printf "%20.15f%20.15f%20.15f\n", psi*57.2957795130823, theta*57.2957795130823, phi*57.2957795130823;
	}
}' "${FNAME}.mov" > Bar2RotC_Euler_3_1_2_MBDyn.dat

awk -v ISTEP="$ISTEP" -v FSTEP="$FSTEP" -v BAR_2_END_NODE="$BAR_2_END_NODE" '$1==BAR_2_END_NODE {
	if (i == 0) {
		print "# R11 R12 R13 R21 R22 R23 R31 R32 R33";
	}
	if (++i > ISTEP && i <= FSTEP) {
		R11 = $5;
		R21 = $8;
		R31 = $11;
		R12 = $6;
		R22 = $9;
		R32 = $12;
		R13 = $7;
		R23 = $10;
		R33 = $13;

		# psi = atan2(-R31, R33);
		# theta = atan2(R32, -R31*sin(psi) + R33*cos(psi));
		# phi = atan2(-R12, R22);
		
		printf "%20.15f%20.15f%20.15f%20.15f%20.15f%20.15f%20.15f%20.15f%20.15f\n", R11, R12, R13, R21, R22, R23, R31, R32, R33;
	}
}' "${FNAME}.mov" > Bar2RotC_R_MBDyn.dat

awk -v ISTEP="$ISTEP" -v FSTEP="$FSTEP" '$1==23 {
	if (i == 0) {
		print "# Rotation[rad]";
	}
	if (++i > ISTEP && i <= FSTEP) {
		phi = -$19*57.2957795130823;
		while (phi - phi0 < -180) {
			phi += 360;
		}
		while (phi - phi0 > 180) {
			phi -= 360;
		}
		printf "%20.15f\n", phi/57.2957795130823;
		phi0 = phi;
	}
}' "${FNAME}.jnt" > RotationC_MBDyn.dat

