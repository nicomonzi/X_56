% MBDyn (C) is a multibody analysis code. 
% http://www.mbdyn.org
% 
% Copyright (C) 1996-2017
% 
% Pierangelo Masarati	<masarati@aero.polimi.it>
% Paolo Mantegazza	<mantegazza@aero.polimi.it>
% 
% Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
% via La Masa, 34 - 20156 Milano, Italy
% http://www.aero.polimi.it
% 
% Changing this copyright notice is forbidden.
% 
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation (version 2 of the License).
% 
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
%
% H.H. Yoo, J. Chung
% "Dynamics of Rectangular Plates Undergoing Prescribed Overall Motion"
% Journal of Sound and Vibration, 239(1), January 4, 2001 (123-137)
% doi:10.1006/jsvi.2000.3111
% 
% Liu Jinyang, Hong Jiazhen
% "Geometric nonlinear formulation and discretization method
% for a rectangular plate undergoing large overall motions"
% Mechanics Research Communications, Volume 32, Issue 5,
% September-October 2005, Pages 561-571
% doi:10.1016/j.mechrescom.2004.10.007
%
% Beerinder Singh, Inderjit Chopra
% "Insect-Based Hover-Capable Flapping Wings for Micro Air Vehicles:
% Experiments and Analysis"
% AIAA Journal, 2008, 46(9), 2115-2135
% doi:10.2514/1.28192

function jinyang2005_5

fd = header('jinyang2005_5_m.set');

%NL = 16;
%NB = 8;
NL = 32;
NB = 16;
%NL = 64;
%NB = 32;

fprintf(fd, 'set: const integer NL = %d;\n', NL);
fprintf(fd, 'set: const integer NB = %d;\n', NB);
fprintf(fd, '\n');
fprintf(fd, 'set: integer ROOT = 0;\n'),
fprintf(fd, 'set: integer PRESCRIBED_OMEGA = 91;\n'),
fprintf(fd, '\n');

fprintf(fd, 'set: const integer NNODES = %d;\n', (NL + 1)*(NB + 1));
fprintf(fd, 'set: const integer NBODIES = %d;\n', (NL + 1)*(NB + 1));
fprintf(fd, 'set: const integer NSHELLS = %d;\n', NL*NB);
fprintf(fd, 'set: const integer NJOINTS = %d;\n', NB + 1);
fprintf(fd, '\n');

E = 7.e+10; % N/m^2
NU = 0.3; % adim
RHO = 3.e+3; % kg/m^3

L = 1.; % m
B = .5; % m
H = 2.5e-3; % m

fprintf(fd, 'set: const real E = %e; # N/m^2\n', E);
fprintf(fd, 'set: const real NU = %e; # adim\n', NU);
fprintf(fd, 'set: const real RHO = %e; # kg/m^3\n', RHO);
fprintf(fd, '\n');

fprintf(fd, 'set: const real L = %e; # m\n', L);
fprintf(fd, 'set: const real B = %e; # m\n', B);
fprintf(fd, 'set: const real H = %e; # m\n', H);
fprintf(fd, '\n');

footer(fd);

fd = header('jinyang2005_5_m.ref');

fprintf(fd, 'reference: ROOT,\n');
fprintf(fd, '        reference, global, null,\n');
fprintf(fd, '        reference, global, eye,\n');
fprintf(fd, '        reference, global, null,\n');
fprintf(fd, '        reference, global, 0., 0., OMEGA_0;\n');

footer(fd);

fd = header('jinyang2005_5_m.nod');

OFFSET = 100;

for il = 0:NL,
	for ib = 0:NB,
		fprintf(fd, 'structural: %d, dynamic,\n', OFFSET*il + ib);
		fprintf(fd, '        reference, ROOT, %e*L, 0., %e*B,\n', il/NL, ib/NB);
		fprintf(fd, '        reference, ROOT, eye,\n');
		fprintf(fd, '        reference, ROOT, null,\n');
		fprintf(fd, '        reference, ROOT, null;\n');
		fprintf(fd, '\n');
	end
end

footer(fd);

fd = header('jinyang2005_5_m.elm');

fprintf(fd, '# joints\n');
fprintf(fd, '\n');

for ib = 0:NB,
	% bug in total pin joint?
	%
	% fprintf(fd, 'joint: %d, total pin joint,\n', ib);
	% fprintf(fd, '        %d,\n', ib);
	% fprintf(fd, '                position, reference, ROOT, null,\n');
	% fprintf(fd, '                position orientation, reference, ROOT, eye,\n');
	% fprintf(fd, '                rotation orientation, reference, ROOT, eye,\n');
	% fprintf(fd, '        position, reference, ROOT, null,\n');
	% fprintf(fd, '        position orientation, reference, ROOT, eye,\n');
	% fprintf(fd, '        rotation orientation, reference, ROOT, eye,\n');
	% fprintf(fd, '        position constraint, 1, 1, 1, null,\n');
	% fprintf(fd, '        orientation constraint, 1, 1, angular velocity, 0., 0., 1., reference, PRESCRIBED_OMEGA;\n');
	% fprintf(fd, '\n');
	fprintf(fd, 'joint: %d, total joint,\n', ib);
	fprintf(fd, '        %d,\n', ib);
	fprintf(fd, '                position, reference, ROOT, null,\n');
	fprintf(fd, '                position orientation, reference, ROOT, eye,\n');
	fprintf(fd, '                rotation orientation, reference, ROOT, eye,\n');
	fprintf(fd, '        9999,\n');
	fprintf(fd, '                position, reference, ROOT, null,\n');
	fprintf(fd, '                position orientation, reference, ROOT, eye,\n');
	fprintf(fd, '                rotation orientation, reference, ROOT, eye,\n');
	fprintf(fd, '        position constraint, 1, 1, 1, null,\n');
	fprintf(fd, '        orientation constraint, 1, 1, angular velocity, 0., 0., 1., reference, PRESCRIBED_OMEGA;\n');
	fprintf(fd, '\n');
end

fprintf(fd, '# rigid bodies\n');
fprintf(fd, '\n');

for il = 0:NL,
	for ib = 0:NB,
		db = (1/NB);
		dl = (1/NL);
		xb = 0;
		xl = 0;
		if (ib == 0),
			db = db/2;
			xb = ib/NB/4;
		elseif (ib == NB),
			db = db/2;
			xb = -ib/NB/4;
		end
		if (il == 0),
			dl = dl/2;
			xl = il/NL/4;
		elseif (il == NL),
			dl = dl/2;
			xl = -il/NL/4;
		end
		m = db*dl;
		Jx = db^2*m/12;
		Jz = dl^2*m/12;

		fprintf(fd, 'body: %d, %d,\n', OFFSET*il + ib, OFFSET*il + ib);
		fprintf(fd, '        %e*B*L*H*RHO,\n', m);
		fprintf(fd, '        reference, node, %e*L, 0., %e*B,\n', xl, xb);
		fprintf(fd, '        diag, %e*B^3*L*H*RHO, (%e*B^3*L + %e*B*L^3)*H*RHO, %e*B*L^3*H*RHO;\n', Jx, Jx, Jz, Jz);
		fprintf(fd, '\n');
	end
end

fprintf(fd, 'inertia: 0, body, all;\n');
fprintf(fd, '\n');

fprintf(fd, '# shells\n');
fprintf(fd, '\n');

for il = 1:NL,
	for ib = 1:NB,
		label = OFFSET*il + ib;
		node1 = OFFSET*il + ib;
		node2 = OFFSET*(il - 1) + ib;
		node3 = OFFSET*(il - 1) + (ib - 1);
		node4 = OFFSET*il + (ib - 1);
		fprintf(fd, 'shell4easans: %d, %d, %d, %d, %d,\n', ...
			label, node1, node2, node3, node4);
		fprintf(fd, '        isotropic, E, E, nu, NU, thickness, H;\n');
	end
end

footer(fd);

endfunction



function fd = header(fname);

fd = fopen(fname, 'w');

fprintf(fd, '# MBDyn (C) is a multibody analysis code. \n');
fprintf(fd, '# http://www.mbdyn.org\n');
fprintf(fd, '# \n');
fprintf(fd, '# Copyright (C) 1996-2017\n');
fprintf(fd, '# \n');
fprintf(fd, '# Pierangelo Masarati	<masarati@aero.polimi.it>\n');
fprintf(fd, '# Paolo Mantegazza	<mantegazza@aero.polimi.it>\n');
fprintf(fd, '# \n');
fprintf(fd, '# Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano\n');
fprintf(fd, '# via La Masa, 34 - 20156 Milano, Italy\n');
fprintf(fd, '# http://www.aero.polimi.it\n');
fprintf(fd, '# \n');
fprintf(fd, '# Changing this copyright notice is forbidden.\n');
fprintf(fd, '# \n');
fprintf(fd, '# This program is free software; you can redistribute it and/or modify\n');
fprintf(fd, '# it under the terms of the GNU General Public License as published by\n');
fprintf(fd, '# the Free Software Foundation (version 2 of the License).\n');
fprintf(fd, '# \n');
fprintf(fd, '# \n');
fprintf(fd, '# This program is distributed in the hope that it will be useful,\n');
fprintf(fd, '# but WITHOUT ANY WARRANTY; without even the implied warranty of\n');
fprintf(fd, '# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n');
fprintf(fd, '# GNU General Public License for more details.\n');
fprintf(fd, '# \n');
fprintf(fd, '# You should have received a copy of the GNU General Public License\n');
fprintf(fd, '# along with this program; if not, write to the Free Software\n');
fprintf(fd, '# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA\n');
fprintf(fd, '#\n');

endfunction



function footer(fd)

fprintf(fd, '\n');
fprintf(fd, '# vim%s\n', ':ft=mbd');

fclose(fd);

endfunction
