# $Header: /var/cvs/mbdyn/mbdyn/mbdyn-1.0/tests/beam/b2.elm,v 1.1 2007/04/15 19:21:05 masarati Exp $

body: 1000+curr_label-1, curr_label-1, 
	dme,
	dle/2, 0., 0.,
	diag, dje, dme/12.*dle^2, dme/12.*dle^2;
body: curr_label, curr_label, 
	dmi,
	null,
	diag, dji, dmi/12*dli^2, dmi/12*dli^2;
body: curr_label+1, curr_label+1, 
	dme,
	-dle/2, 0., 0.,
	diag, dje, dme/12.*dle^2, dme/12.*dle^2;

beam: curr_label,
	curr_label-1, null,
	curr_label,   null,
	curr_label+1, null,
	eye,
	linear elastic generic,
	diag, EA, GAy, GAz, GJ, EJy, EJz,
	same,
	same;
