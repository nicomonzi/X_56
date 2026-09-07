# $Header: /var/cvs/mbdyn/mbdyn/mbdyn-1.0/tests/beam/DETC2007-35511/beam.elm,v 1.2 2009/08/21 12:04:29 masarati Exp $
#
# Copyright 2009 Pierangelo Masarati <masarati@aero.polimi.it>
#
# Space manipulator based on DETC2007-35511
# Proceedings of IDETC/CIE 2007
# ASME 2007 International Design Engineering Technical Conferences &
# Computers and Information in Engineering Conference
# September 4-7, 2007, Las Vegas, USA
# INVESTIGATION OF BOUNDARY CONDITIONS FOR FLEXIBLE MULTIBODY
# SPACECRAFT DYNAMICS (DRAFT VERSION)
# John R. MacLean
# An Huynh
# Leslie J. Quiocho

body: CURR_BEAM + 1000, PIN + 2*CURR_BEAM - 2,
	DM/2.,
	reference, node, DL/4., 0., 0.,
	diag,
		DJ/2., (DL/2.)^2*DM/2./12., (DL/2.)^2*DM/2./12.;

body: CURR_BEAM + 2000, PIN + 2*CURR_BEAM - 1,
	DM,
	reference, node, null,
	diag,
		DJ, DL^2*DM/12., DL^2*DM/12.;

body: CURR_BEAM + 3000, PIN + 2*CURR_BEAM - 0,
	DM/2.,
	reference, node, -DL/4., 0., 0.,
	diag,
		DJ/2., (DL/2.)^2*DM/2./12., (DL/2.)^2*DM/2./12.;

beam3: CURR_BEAM,
	PIN + 2*CURR_BEAM - 2, null,
	PIN + 2*CURR_BEAM - 1, null,
	PIN + 2*CURR_BEAM - 0, null,
	reference, PIN, eye,
	linear elastic generic,
	# linear viscoelastic generic,
		diag, EA, GA, GA, GJ, EIY, EI,
		# proportional, 1e-5,
	same,
	same;

