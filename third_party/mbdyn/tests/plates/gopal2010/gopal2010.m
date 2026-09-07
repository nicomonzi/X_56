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
% Based on Pradeep Gopalakrishnan and Danesh K. Tafti
% "Effect of wing flexibility on lift and thrust production in flapping flight."
% AIAA Journal, 48(5):865-877, May 2010. doi:10.2514/1.39957

function gopal2010

fd = header('gopal2010_m.set');

% Ria's mail 10-14-2011
% 
% chord = 2cm,
% span (including root cut out) = 4*2 = 8cm
% thickness = 0.05 mm
% U_infinity = 4m/s
% f = 47.75 Hz
% rho_polyester =  1350 kg/m3
%
% Grid motion:
% alpha = pitch_knot*pi/180+pitch_angle*pi/180*cos(omega_t+phase*pi/180)
% flap_angle = flap_amp*pi/180*(cos(omega_t))
% pitch_knot = 12.5 deg
% pitch_angle = 32.5 deg
% phase = 90 deg
% flap_amp = 30 deg
%
% Ria's mail 10-17-2011
%
% I think while making the grid, I have non dimensionalized it wrt c=4cm,
% which is why I have t/c=0.00125.
%
% Sorry about that, could you please use these parameters,
% C=4cm,
% f = 11.94 Hz
% U_infty = 2m/s

E = 2e+9; % N/m^2; polyethilene; 3.5 GPa from the internet
NU = 0.25; % adim; polyethilene, from the internet
RHO = 1350.; % kg/m^3

PRESTRESS = 0.; % N/m

B = 40e-3; % m, chord
L = 3.5*B; % m
H = .05e-3; % m
L0 = 0.5*B; % m

FREQ = 11.94; % Hz
PITCH_0 = 12.5; % deg
PITCH_A = 32.5; % deg
FLAP_A = 30.; % deg
PHASE = 90.; % deg

NB = 10;
% NL = 3.5*NB;
NL = 2*NB;

fprintf(fd, 'set: const integer NL = %d;\n', NL);
fprintf(fd, 'set: const integer NB = %d;\n', NB);
fprintf(fd, '\n');
fprintf(fd, 'set: integer GROUND = 99997;\n'),
fprintf(fd, 'set: integer FLAP = 99998;\n'),
fprintf(fd, 'set: integer PITCH = 99999;\n'),
fprintf(fd, '\n');

fprintf(fd, 'set: const integer NNODES = %d;\n', (NL + 1)*(NB + 1));
fprintf(fd, 'set: const integer NBODIES = %d;\n', (NL + 1)*(NB + 1));
fprintf(fd, 'set: const integer NSHELLS = %d;\n', NL*NB);
fprintf(fd, 'set: const integer NJOINTS = %d;\n', 2*(NL + NB));
fprintf(fd, '\n');

fprintf(fd, 'set: real E = %e; # N/m^2; tune as needed\n', E);
fprintf(fd, 'set: const real NU = %e; # adim\n', NU);
fprintf(fd, 'set: real RHO = %e; # kg/m^3; tune as needed\n', RHO);
fprintf(fd, '\n');

fprintf(fd, 'set: const real L0 = %e; # m\n', L0);
fprintf(fd, 'set: const real L = %e; # m\n', L);
fprintf(fd, 'set: const real B = %e; # m\n', B);
fprintf(fd, 'set: real H = %e; # m; tune as needed\n', H);
fprintf(fd, 'set: real PS = %e; # m\n', PRESTRESS);
fprintf(fd, '\n');

fprintf(fd, 'set: real FREQ = %e*2*pi; # Hz\n', FREQ);
fprintf(fd, 'set: real PITCH_0 = %e*deg2rad; # deg\n', PITCH_0);
fprintf(fd, 'set: real PITCH_A = %e*deg2rad; # deg\n', PITCH_A);
fprintf(fd, 'set: real FLAP_A = %e*deg2rad; # deg\n', FLAP_A);
fprintf(fd, 'set: real PHASE = %e*deg2rad; # deg\n', PHASE);
fprintf(fd, '\n');

footer(fd);

fd = header('gopal2010_m.ref');

fprintf(fd, 'reference: GROUND,\n');
fprintf(fd, '        reference, global, null,\n');
fprintf(fd, '        reference, global, eye,\n');
fprintf(fd, '        reference, global, null,\n');
fprintf(fd, '        reference, global, null;\n');

fprintf(fd, 'reference: FLAP,\n');
fprintf(fd, '        reference, GROUND, null,\n');
fprintf(fd, '        reference, GROUND, eye,\n');
fprintf(fd, '        reference, GROUND, null,\n');
fprintf(fd, '        reference, GROUND, null;\n');

fprintf(fd, 'reference: PITCH,\n');
fprintf(fd, '        reference, FLAP, 0., L0, 0.,\n');
fprintf(fd, '        reference, FLAP, eye,\n');
fprintf(fd, '        reference, FLAP, null,\n');
fprintf(fd, '        reference, FLAP, null;\n');

footer(fd);

fd = header('gopal2010_m.nod');

OFFSET = 100;

for il = 0:NL,
	for ib = 0:NB,
		fprintf(fd, 'structural: %d, dynamic,\n', OFFSET*il + ib);
		fprintf(fd, '        reference, PITCH, %e*B, %e*L, 0.,\n', ib/NB - 0.25, il/NL);
		fprintf(fd, '        reference, PITCH, eye,\n');
		fprintf(fd, '        reference, PITCH, null,\n');
		fprintf(fd, '        reference, PITCH, null;\n');
		fprintf(fd, '\n');
	end
end

footer(fd);

fd = header('gopal2010_m.elm');

fprintf(fd, '# joints\n');
fprintf(fd, '\n');

boundary_nodes = [[1:NL]*OFFSET, [1:NB] + NL*OFFSET, [NL - 1:-1:0]*OFFSET + NB, [NB - 1:-1:0]];

for ib = 1:length(boundary_nodes),
	fprintf(fd, 'joint: %d, total joint,\n', boundary_nodes(ib));
	fprintf(fd, '        PITCH,\n');
	fprintf(fd, '                position, reference, PITCH, null,\n');
	fprintf(fd, '                position orientation, reference, PITCH, eye,\n');
	fprintf(fd, '                rotation orientation, reference, PITCH, eye,\n');
	fprintf(fd, '        %d,\n', boundary_nodes(ib));
	fprintf(fd, '                position, reference, PITCH, null,\n');
	fprintf(fd, '                position orientation, reference, PITCH, eye,\n');
	fprintf(fd, '                rotation orientation, reference, PITCH, eye,\n');
	fprintf(fd, '        position constraint, 1, 1, 1, null,\n');
	fprintf(fd, '        orientation constraint, 1, 1, 1, null;\n');
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
		fprintf(fd, '        reference, node, %e*B, %e*L, 0.,\n', xb, xl);
		fprintf(fd, '        diag, %e*B*L^3*H*RHO, %e*B^3*L*H*RHO, (%e*B^3*L + %e*B*L^3)*H*RHO;\n', Jz, Jx, Jz, Jx);
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
		fprintf(fd, '        isotropic, E, E, nu, NU, thickness, H,\n');
		fprintf(fd, '        prestress, PS, 0., 0., 0., PS, 0., 0., 0., 0., 0., 0., 0.;\n');
	end
end

footer(fd);

fd = header('gopal2010_ext_m.elm');

mb_file = 'gopal2010_shell_mb.dat';
cfd_file = 'gopal2010_shell_surface_grid.dat';
mapping_file = 'gopal2010_shell_H.dat';

fprintf(fd, 'force: 111, external structural mapping,\n');
fprintf(fd, '        socket,\n');
fprintf(fd, '        create, yes,\n');
fprintf(fd, '        path, "$MBSOCK",\n');
fprintf(fd, '        no signal,\n');
fprintf(fd, '        coupling, tight,\n');
fprintf(fd, '        reference node, PITCH,\n');
fprintf(fd, '        orientation, orientation matrix,\n');
fprintf(fd, '        use reference node forces, no,\n');
fprintf(fd, '        points number, 2*NNODES,\n');
for il = 0:NL,
	for ib = 0:NB,
		fprintf(fd, '                %4d, offset, 0., 0., .25*H, offset, 0., 0., -.25*H,\n', OFFSET*il + ib);
	end
end
fprintf(fd, '        # echo, "%s", precision, 16, surface, "%s", output, "%s", order, 2, basenode, 12, weight, 2, stop;\n', mb_file, cfd_file, mapping_file);
fprintf(fd, '        mapped points number, from file,\n');
fprintf(fd, '        sparse mapping file, "%s";\n', mapping_file);

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
