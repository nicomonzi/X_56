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
}

BEGIN {
	nodefile = "cook_membrane.nod";
	elemfile = "cook_membrane.elm";

	#
	#   1 +
	#     |\
	# 2bs | \
	#     |  \
	#     |   \
	#  11 +    + 111
	#      \   | 2bt
        #       \  |
        #     |  \ |
	#     |   \|
	#     |    + 121
	#     |  s |
	#     +----+

	s  = 48;
	bs = 44;
	bt = 16;

	x1 = 0.;
	x11 = x1 + bs;
	x111 = x1 + bs;
	x121 = x111 + bt;
	y1 = 0.;
	y11 = 0.;
	y111 = s;
	y121 = s;

	th = 1e0;

	E  = 1;
	nu = 0.333;

	if (NCHORD == 0) {
		NCHORD = 8;
	}

	if (NSPAN == 0) {
		NSPAN = 8;
	}

	# NODES
	header(nodefile);
	for (ir = 1; ir <= (NSPAN + 1); ir++) {
		y = y1 + (y111 - y1)*(ir - 1)/NSPAN;
		for (ic = 1; ic <= (NCHORD + 1); ic++) {
			type = "static";
			# if (ir == 1) {
			# 	type = "static";
			# } else {
			# 	type = "static"; # displacement";
			# }
			print "structural: " ((NCHORD + 1)*(ir - 1) + ic) ", " type "," >> nodefile;
			xr1 = x1 + (x11 - x1)*(ic - 1)/NCHORD;
			xr11 = x111 + (x121 - x111)*(ic - 1)/NCHORD;
			x = xr1 + (xr11 - xr1)*(ir - 1)/NCHORD;
			printf "\t%.16e, %.16e, %.16e,\n", x, y, 0. >> nodefile;
			if (type == "static") {
				print "\teye," >> nodefile;
				print "\tnull," >> nodefile;
			}
			print "\tnull;" >> nodefile;
			print "" >> nodefile;
		}
	}

	printf "# %sft=mbd\n", "vim:" >> nodefile;

	# elements
	header(elemfile);
	print "# joints" >> elemfile;
	for (ir = 1; ir <= (NSPAN + 1); ir++) {
		for (ic = 1; ic <= (NCHORD + 1); ic++) {
			printf "joint: %d, total pin joint,\n", (NCHORD + 1)*(ir - 1) + ic >> elemfile;
			printf "\t%d,\n", (NCHORD + 1)*(ir - 1) + ic >> elemfile;
			print "\t\tposition, reference, global, null," >> elemfile;
			print "\t\tposition orientation, reference, global, eye," >> elemfile;
			print "\t\tposition, reference, global, null," >> elemfile;
			print "\t\tposition orientation, reference, global, eye," >> elemfile;
			if ( ir == 1) {
				print "\tposition constraint, active, active, active," >> elemfile;
			} else {
				print "\tposition constraint, inactive, inactive, active," >> elemfile;
			}
			print "\t\tnull," >> elemfile;
			print "\torientation constraint, active, active, active," >> elemfile;
			print "\t\tnull;" >> elemfile;
		}
	}


	# shell = "shell4easans";
	shell = "membrane4eas";
	print "# shell elements" >> elemfile;
	for (ir = 1; ir <= NSPAN; ir++) {
		chord = (x11 - x1)*(NSPAN + 1 - ir)/NSPAN + (x121 - x111)*(ir - 1)/NSPAN;
		for (ic = 1; ic <= NCHORD; ic++) {
			print shell ": " (NCHORD*(ir - 1) + ic) "," >> elemfile;
			print "\t" ((NCHORD + 1)*ir + ic) ", " ((NCHORD + 1)*(ir - 1) + ic) ", " ((NCHORD + 1)*(ir - 1) + ic + 1) ", " ((NCHORD + 1)*ir + ic + 1) "," >> elemfile;
			printf "\tisotropic, E, %.16e, nu, %.16e, thickness, %.16e;\n", E, nu, th >> elemfile;
			print "" >> elemfile;
		}
	}

	ir = (NSPAN + 1)
		for (ic = 1; ic <= (NCHORD + 1); ic++) {
			print "force: " ((NCHORD + 1)*(ir - 1) + ic) ", absolute displacement," >> elemfile;
			print "\t" ((NCHORD + 1)*(ir - 1) + ic) ", " >> elemfile;
			# print "\tposition, null," >> elemfile;
			print "\t1., 0., 0.," >> elemfile;
			printf "\t\tcosine, 0., 2*pi/1, %.8e, half, 0.;", 0.5 * 1e0/(NCHORD+1) >> elemfile;
			print "" >> elemfile;
		}

	printf "# %sft=mbd", "vim:" >> elemfile;
}

