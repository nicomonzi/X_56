# $Header: /var/cvs/mbdyn/mbdyn/mbdyn-1.0/tests/beam/b.elm,v 1.5 2007/04/15 16:45:11 masarati Exp $

body: curr_label, curr_label, 
      dm,
      null,
      diag, dj, dm/12*(dl)^2, dm/12*(dl)^2;
body: curr_label+1, curr_label+1, 
      dm,
      null,
      diag, dj, dm/12*(dl)^2, dm/12*(dl)^2;

beam: curr_label,
      curr_label-1, null,
      curr_label,   null,
      curr_label+1, null,
      eye,
      linear elastic generic,
      diag, EA, GAy, GAz, GJ, EJy, EJz,
      same,
      same;
