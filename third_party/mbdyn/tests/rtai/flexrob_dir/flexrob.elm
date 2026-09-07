# $Header: /var/cvs/mbdyn/mbdyn/mbdyn-1.0/tests/rtai/flexrob_dir/flexrob.elm,v 1.3 2007/10/03 14:02:23 masarati Exp $
# Prepared by Michele Attolico

body: curr_label, curr_label, 
      Dm,
      null,
      diag, jxx, jyy, jzz;

beam: curr_label,
      	curr_label-10, null,
      	curr_label,   null,
      	curr_label+10, null,
	reference, body, eye,
	linear viscoelastic generic,
		diag,
		E*A,
		G*Ay,
		G*Az,
		G*Jx,
		E*Jy,
		E*Jz,
		proportional, beam_damp,
		same,
		same;
