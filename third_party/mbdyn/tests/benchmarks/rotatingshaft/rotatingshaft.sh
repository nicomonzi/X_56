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

# run the analysis as
#	$ octave -q rotatingshaft_gen.m
#	$ mbdyn -f rotatingshaft -o x
# with 16 beam elements and a time step of 1e-3, starting from -10 s
#
# run 
#	$ awk -f abs2rel.awk -v InputMode=matrix -v RefNode=0 x.mov > xr.mov

OUTFNAME=x
OUTFNAMEREL=xr

# Shaft mid-point displacement components
awk '$1 == 16 && i++ >= 10000 {print (1e-3*(i-1) - 10), $3}' ${OUTFNAME}.mov > MBDyn_BeamMidDisplacements_u2.dat
awk '$1 == 16 && i++ >= 10000 {print (1e-3*(i-1) - 10), $4}' ${OUTFNAME}.mov > MBDyn_BeamMidDisplacements_u3.dat

# Shaft mid-point velocity components
awk '$1 == 16 && i++ >= 10000 {print (1e-3*(i-1) - 10), $15}' ${OUTFNAME}.mov > MBDyn_BeamMidVelocities_v2.dat
awk '$1 == 16 && i++ >= 10000 {print (1e-3*(i-1) - 10), $16}' ${OUTFNAME}.mov > MBDyn_BeamMidVelocities_v3.dat

# Shaft mid-point shear force components
#awk '$1 == 8 && i++ >= 10000 {f2 = $9; getline; f2 += $3; print (1e-3*(i-1) - 10), f2/2}' ${OUTFNAME}.act > MBDyn_BeamMidForces_F2.dat
#awk '$1 == 8 && i++ >= 10000 {f3 = $10; getline; f3 += $4; print (1e-3*(i-1) - 10), f3/2}' ${OUTFNAME}.act > MBDyn_BeamMidForces_F3.dat
awk '$1 == 8 && i++ >= 10000 {print (1e-3*(i-1) - 10), ((1 - sqrt(3))*$3 + (1 + sqrt(3))*$9)/2}' ${OUTFNAME}.act > MBDyn_BeamMidForces_F2.dat
awk '$1 == 8 && i++ >= 10000 {print (1e-3*(i-1) - 10), ((1 - sqrt(3))*$4 + (1 + sqrt(3))*$10)/2}' ${OUTFNAME}.act > MBDyn_BeamMidForces_F3.dat

# Shaft mid-point bending moment components
awk '$1 == 8 && i++ >= 10000 {print (1e-3*(i-1) - 10), ((1 - sqrt(3))*$6 + (1 + sqrt(3))*$12)/2}' ${OUTFNAME}.act > MBDyn_BeamMidMoments_M2.dat
awk '$1 == 8 && i++ >= 10000 {print (1e-3*(i-1) - 10), ((1 - sqrt(3))*$7 + (1 + sqrt(3))*$13)/2}' ${OUTFNAME}.act > MBDyn_BeamMidMoments_M3.dat

# Trajectory of mid-point
awk '$1 == 16 && i++ >= 10000 {print $3, $4}' ${OUTFNAME}.mov > MBDyn_BeamMidTrajectoryAbs.dat
awk '$1 == 16 && i++ >= 10000 {print $3, $4}' ${OUTFNAMEREL}.mov > MBDyn_BeamMidTrajectoryRel.dat

# Phase plane of mid-point
awk '$1 == 16 && i++ >= 10000 {print $3, $15}' ${OUTFNAME}.mov > MBDyn_BeamMidPhasePlane_u2v2.dat
awk '$1 == 16 && i++ >= 10000 {print $4, $16}' ${OUTFNAME}.mov > MBDyn_BeamMidPhasePlane_u3v3.dat
