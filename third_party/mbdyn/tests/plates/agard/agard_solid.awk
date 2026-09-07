# $Header: /var/cvs/mbdyn/mbdyn/mbdyn-1.0/tests/plates/agard/agard_solid.awk,v 1.10 2017/01/12 15:06:33 masarati Exp $

function header(out) {
	printf "" > out;
	print "# MBDyn (C) is a multibody analysis code. " >> out;
	print "# http://www.mbdyn.org" >> out;
	print "# " >> out;
	print "# Copyright (C) 1996-2017" >> out;
	print "# " >> out;
	print "# Pierangelo Masarati	<masarati@aero.polimi.it>" >> out;
	print "# Paolo Mantegazza	<mantegazza@aero.polimi.it>" >> out;
	print "# " >> out;
	print "# Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano" >> out;
	print "# via La Masa, 34 - 20156 Milano, Italy" >> out;
	print "# http://www.aero.polimi.it" >> out;
	print "# " >> out;
	print "# Changing this copyright notice is forbidden." >> out;
	print "# " >> out;
	print "# This program is free software; you can redistribute it and/or modify" >> out;
	print "# it under the terms of the GNU General Public License as published by" >> out;
	print "# the Free Software Foundation (version 2 of the License)." >> out;
	print "# " >> out;
	print "# " >> out;
	print "# This program is distributed in the hope that it will be useful," >> out;
	print "# but WITHOUT ANY WARRANTY; without even the implied warranty of" >> out;
	print "# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the" >> out;
	print "# GNU General Public License for more details." >> out;
	print "# " >> out;
	print "# You should have received a copy of the GNU General Public License" >> out;
	print "# along with this program; if not, write to the Free Software" >> out;
	print "# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA" >> out;
	print "# " >> out;
	print "# AGARD 445.6 wing model" >> out;
	print "# Author: Pierangelo Masarati" >> out;
	print "#         based on the Code Aster model provided by Giulio Romanelli" >> out;
}

BEGIN {
	nodefile = "agard_solid.dat";

	#   0.000000000000000   0.000000000000000
	#   0.557784000000000   0.000000000000000
	#   0.809396000000000   0.762000000000000
	#   1.177595000000000   0.762000000000000
	#x1   = 0;
	#x11  = 0.557784000000000;
	#x111 = 0.809396000000000;
	#x121 = 1.177595000000000;
	#y1   = 0;
	#y11  = 0;
	#y111 = 0.762000000000000;
	#y121 = 0.762000000000000;

	s = 2.5 * 0.3048;	# ft -> m
	bs = 1.833 * 0.3048;	# ft -> m
	bt = 1.208 * 0.3048;	# ft -> m

	x1 = 0.;
	x11 = x1 + bs;
	x111 = x1 + bs/4. + s - bt/4.;
	x121 = x111 + bt;
	y1 = 0.;
	y11 = 0.;
	y111 = s;
	y121 = s;

	xp[1] = 0.0;
	xp[2] = .005;
	xp[3] = .0075;
	xp[4] = .0125;
	xp[5] = .025;
	xp[6] = .050;
	xp[7] = .075;
	xp[8] = .10;
	xp[9] = .15;
	xp[10] = .20;
	xp[11] = .25;
	xp[12] = .30;
	xp[13] = .35;
	xp[14] = .40;
	xp[15] = .45;
	xp[16] = .50;
	xp[17] = .55;
	xp[18] = .60;
	xp[19] = .65;
	xp[20] = .70;
	xp[21] = .75;
	xp[22] = .80;
	xp[23] = .85;
	xp[24] = .90;
	xp[25] = .95;
	xp[26] = 1.00;

	zp[1] = 0.00;
	zp[2] = .00311;
	zp[3] = .00378;
	zp[4] = .00481;
	zp[5] = .00656;
	zp[6] = .00877;
	zp[7] = .01062;
	zp[8] = .01216;
	zp[9] = .01463;
	zp[10] = .01649;
	zp[11] = .01790;
	zp[12] = .01894;
	zp[13] = .01962;
	zp[14] = .01996;
	zp[15] = .01996;
	zp[16] = .01952;
	zp[17] = .01867;
	zp[18] = .01742;
	zp[19] = .01584;
	zp[20] = .01400;
	zp[21] = .01193;
	zp[22] = .00966;
	zp[23] = .00728;
	zp[24] = .00490;
	zp[25] = .00249;
	zp[26] = 0.0; # .00009

	NCHORD = 25;
	NSPAN = 20;

	NPOINTS = (NSPAN + 1)*2*NCHORD;

	header(nodefile);
	printf "%d\n", NPOINTS >> nodefile;
	for (ir = 1; ir <= (NSPAN + 1); ir++) {
		y = y1 + s*(ir - 1)/NSPAN;
		for (ic = 1; ic <= (NCHORD + 1); ic++) {
			xle = x1 + (x111 - x1)*(ir - 1)/NSPAN;
			chord = bs*(NSPAN + 1 - ir)/NSPAN + bt*(ir - 1)/NSPAN;
			x = xle + xp[ic]*chord;
			printf "%+15.7e%+15.7e%+15.7e\n", x, y, zp[ic]*chord >> nodefile;
			if (ic > 1 && ic <= NCHORD) {
				printf "%+15.7e%+15.7e%+15.7e\n", x, y, -zp[ic]*chord >> nodefile;
			}
		}
	}
}
