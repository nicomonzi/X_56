# $Header: /var/cvs/mbdyn/mbdyn/mbdyn-1.0/tests/plates/agard/agard.awk,v 1.20 2017/01/12 15:06:33 masarati Exp $

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
	print "" >> out;
}

BEGIN {
	nodefile = "agard_m.nod";
	elemfile = "agard_m.elm";

	# From NASA-TM-100492
	# <http://hdl.handle.net/2060/19880001820>

	# Figure 1. - Case #2, 
	# Structure	Mounting	s		2bs, ft		2bt, ft		Lambda, deg
	# Solid		Wall		2.5		1.833 (11/6)	1.208		43.15		ft
	#				0.76200		0.55870		0.36820				m

	#   x, chordwise, m	y, spanwise, m
	#   0.000000000000000   0.000000000000000
	#   0.557784000000000   0.000000000000000
	#   0.809396000000000   0.762000000000000
	#   1.177595000000000   0.762000000000000
	#
	#   1 +
	#     |\
	# 2bs | \
	#     |  \
	#     |   + 111
	#  11 +   |
	#      \  | 2bt
        #       \ |
        #     |  \|
	#     |   + 121
	#     | s |
	#     +---+

	# x1   = 0;
	# x11  = 0.557784000000000;
	# x111 = 0.809396000000000;
	# x121 = 1.177595000000000;
	# y1   = 0;
	# y11  = 0;
	# y111 = 0.762000000000000;
	# y121 = 0.762000000000000;

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

	z[1] = 0.000000;
	z[2] = 0.012160;
	z[3] = 0.016490;
	z[4] = 0.018940;
	z[5] = 0.019960;
	z[6] = 0.019520;
	z[7] = 0.017420;
	z[8] = 0.014000;
	z[9] = 0.009660;
	z[10] = 0.004900;
	z[11] = 0.000090;

	# Table I of TN D-1616 (appendix of AGARD report)
	# span 2.5 ft, "Weakened", Model 2: 0.13665 slugs
	# span 2.5 ft, "Weakened", Model 5: 0.12143 slugs
	total_mass = 0.12143 * 14.5939029372064;	# slugs -> kg
	airfoil_surface = (z[1] + z[11])/2;
	for (i = 2; i < 11; i++) {
		airfoil_surface += z[i];
	}
	airfoil_surface *= 2./10.;
	volume = airfoil_surface*(bs + bt)/2*(bs + bt)/2*s;
	rho = total_mass/volume;
	### solid
	# E = 1.1846e6 * .45359237*9.80665/(.0254*.0254);	# psi -> N/m^2
	# G = 0.1159e6 * .45359237*9.80665/(.0254*.0254);	# psi -> N/m^2
	### weakened
	# E = 0.47072e6 * .45359237*9.80665/(.0254*.0254);	# psi -> N/m^2
	# G = 0.059745e6 * .45359237*9.80665/(.0254*.0254);	# psi -> N/m^2
	# nu = 0.310;					# non-dim.
	### weakened - Giulio Romanelli
	E_l = 3.151E9;
	E_t = 4.162E8;
	G = 4.392E8;
	nu_lt = 0.31;

	if (NCHORD == 0) {
		NCHORD = 10;
	}

	if (NSPAN == 0) {
		NSPAN = 10;
	}

	# nodes
	header(nodefile);
	for (ir = 1; ir <= (NSPAN + 1); ir++) {
		y = y1 + (y111 - y1)*(ir - 1)/NSPAN;
		for (ic = 1; ic <= (NCHORD + 1); ic++) {
			if (ir == 1) {
				type = "static";
			} else {
				type = "dynamic";
			}
			print "structural: " ((NCHORD + 1)*(ir - 1) + ic) ", " type "," >> nodefile;
			xr1 = x1 + (x11 - x1)*(ic - 1)/NCHORD;
			xr11 = x111 + (x121 - x111)*(ic - 1)/NCHORD;
			x = xr1 + (xr11 - xr1)*(ir - 1)/NCHORD;
			printf "\t%.16e, %.16e, %.16e,\n", x, y, 0. >> nodefile;
			print "\teye," >> nodefile;
			print "\tnull," >> nodefile;
			print "\tnull;" >> nodefile;
			print "" >> nodefile;
		}
	}

	printf "# %sft=mbd\n", "vim:" >> nodefile;

	# elements
	header(elemfile);
	print "# joints" >> elemfile;
	for (ir = 3; ir <= (NCHORD - 1); ir++) {
		print "joint: " ir ", clamp, " ir ", node, node;" >> elemfile;
		print "" >> elemfile;
	}

	V = 0.;

	print "# rigid bodies" >> elemfile;
	DY = (y111 - y1)/NSPAN;
	YCG = 0.;
	for (ir = 2; ir <= (NSPAN + 1); ir++) {
		if (ir == NSPAN + 1) {
			YCG = -DY/4;
			DY /= 2;
		}

		chord = (x11 - x1)*(NSPAN + 1 - ir)/NSPAN + (x121 - x111)*(ir - 1)/NSPAN;

		DX = chord/NCHORD;

		# print "ir=" ir " chord=" chord " DY=" DY " DX=" DX >> "/dev/stdout";

		for (ic = 1; ic <= (NCHORD + 1); ic++) {
			factor = 1.;
			if (ic == 1) {
				h = 2*chord*(z[ic] + z[ic + 1])/2;
				factor = .5;
			} else if (ic == NCHORD + 1) {
				h = 2*chord*(z[ic-1] + z[ic])/2;
				factor = .5;
			} else {
				h = 2*chord*(z[ic - 1] + 2*z[ic] + z[ic + 1])/4;
			}
			DV = factor*DX*DY*h;
			V += DV;
			DM = DV*rho;
			print "body: " ((NCHORD + 1)*(ir - 1) + ic) ", " ((NCHORD + 1)*(ir - 1) + ic) "," >> elemfile;
			printf "\t%.16e,\n", DM >> elemfile;
			if (ic == 1) {
				XCG = DX/4;
			} else if (ic == NCHORD + 1) {
				XCG = -DX/4;
			} else {
				XCG = 0.;
			}
			if (XCG != 0. || YCG != 0.) {
				printf "\t%.16e, %.16e, %.16e,\n", XCG, YCG, 0.  >> elemfile;
			} else {
				print "\tnull," >> elemfile;
			}
			print "\tdiag, 1e-9, 1e-9, 1e-9;" >> elemfile;
			print "" >> elemfile;
		}
	}

	printf "remark: \"total mass\", %.16e;\n", total_mass >> elemfile;
	print "inertia: 0, body, all;" >> elemfile;
	print "" >> elemfile;

	# print "V=" V "; volume=" volume >> "/dev/stdout";

	shell = "shell4eas";
	print "# shell elements" >> elemfile;
	for (ir = 1; ir <= NSPAN; ir++) {
		chord = (x11 - x1)*(NSPAN + 1 - ir)/NSPAN + (x121 - x111)*(ir - 1)/NSPAN;
		for (ic = 1; ic <= NCHORD; ic++) {
			h = 2*chord*(z[ic] + z[ic + 1])/2;
			print shell ": " (NCHORD*(ir - 1) + ic) "," >> elemfile;
			print "\t" ((NCHORD + 1)*ir + ic) ", " ((NCHORD + 1)*(ir - 1) + ic) ", " ((NCHORD + 1)*(ir - 1) + ic + 1) ", " ((NCHORD + 1)*ir + ic + 1) "," >> elemfile;
			printf "\tplane stress orthotropic, E_l, %.16e, E_t, %.16e, nu_lt, %.16e, G, %.16e, thickness, %.16e, as, 1, at, 0.01;\n", E_l, E_t, nu_lt, G, h >> elemfile;
			print "" >> elemfile;
		}
	}

	printf "# %sft=mbd", "vim:" >> elemfile;
}

