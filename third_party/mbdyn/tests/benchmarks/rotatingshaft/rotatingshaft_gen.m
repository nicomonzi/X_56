% MBDyn (C) is a multibody analysis code. 
% http://www.mbdyn.org
% 
% Copyright (C) 1996-2021
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

N_BEAMS = 16;
N_NODES = 2*N_BEAMS + 1;

L = 6; % m
rho = 7800; % kg/m^3
E = 210e9; % Pa
nu = 0.3;
rI = 0.045; % m
rO = 0.05; % m
m = pi*rho*(rO^2 - rI^2); % 11.64; % kg/m
JXX = pi/4*rho*(rO^4 - rI^4);
J11 = 2*JXX;
% JXX = 0;

fd = fopen('rotatingshaft_m.set', 'w');

fprintf(fd, '# MBDyn (C) is a multibody analysis code.\n');
fprintf(fd, '# http://www.mbdyn.org\n');
fprintf(fd, '# \n');
fprintf(fd, '# Copyright (C) 1996-2021\n');
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
fprintf(fd, '\n');

fprintf(fd, 'set: const integer N_BEAMS = %d;\n', N_BEAMS);
fprintf(fd, 'set: const integer N_NODES = %d;\n', N_NODES);

fprintf(fd, '\n');

fprintf(fd, 'set: const integer SHAFT_R = 0;\n');
fprintf(fd, 'set: const integer SHAFT_T = 2*N_BEAMS;\n');

fprintf(fd, '\n');

fprintf(fd, 'set: const real L = %e;\n', L);
fprintf(fd, 'set: const real EA = 313.4e6;\n');
fprintf(fd, 'set: const real GA = 60.5e6;\n');
fprintf(fd, 'set: const real GJ = 272.7e3;\n');
fprintf(fd, 'set: const real EJ = 354.5e3;\n');
fprintf(fd, 'set: real DAMPING_FACTOR = 0.;\n');

fprintf(fd, '\n');

fprintf(fd, 'set: const real A1 = 0.8;\n');
fprintf(fd, 'set: const real A2 = 1.2;\n');
fprintf(fd, 'set: const real T0 = 0.;\n');
fprintf(fd, 'set: const real T1 = 0.5;\n');
fprintf(fd, 'set: const real T2 = 1.;\n');
fprintf(fd, 'set: const real T3 = 1.25;\n');
fprintf(fd, 'set: const real OMEGA = 60.;\n');

fclose(fd);



fd = fopen('rotatingshaft_m.ref', 'w');

fprintf(fd, '# MBDyn (C) is a multibody analysis code.\n');
fprintf(fd, '# http://www.mbdyn.org\n');
fprintf(fd, '# \n');
fprintf(fd, '# Copyright (C) 1996-2021\n');
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
fprintf(fd, '\n');

fprintf(fd, 'reference: SHAFT_R,\n');
fprintf(fd, '\treference, global, null,\n');
fprintf(fd, '\treference, global, eye,\n');
fprintf(fd, '\treference, global, null,\n');
fprintf(fd, '\treference, global, null;\n');

fprintf(fd, 'reference: SHAFT_T,\n');
fprintf(fd, '\treference, global, L, 0., 0.,\n');
fprintf(fd, '\treference, global, eye,\n');
fprintf(fd, '\treference, global, null,\n');
fprintf(fd, '\treference, global, null;\n');

fprintf(fd, '\n');

fprintf(fd, 'constitutive law: 6, 6, linear viscoelastic generic,\n');
fprintf(fd, '\tdiag, EA, GA, GA, GJ, EJ, EJ, proportional, DAMPING_FACTOR;\n');

fclose(fd);



fd = fopen('rotatingshaft_m.nod', 'w');

fprintf(fd, '# MBDyn (C) is a multibody analysis code.\n');
fprintf(fd, '# http://www.mbdyn.org\n');
fprintf(fd, '# \n');
fprintf(fd, '# Copyright (C) 1996-2021\n');
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
fprintf(fd, '\n');

for i = 0:2*N_BEAMS,
	fprintf(fd, 'structural: %d, dynamic,\n', i);
	fprintf(fd, '\treference, SHAFT_R, %d*L/(2*N_BEAMS), 0., 0.,\n', i);
	fprintf(fd, '\treference, SHAFT_R, eye,\n');
	fprintf(fd, '\treference, SHAFT_R, null,\n');
	fprintf(fd, '\treference, SHAFT_R, null;\n');

	fprintf(fd, '\n');
end

fclose(fd);



fd = fopen('rotatingshaft_m.elm', 'w');

fprintf(fd, '# MBDyn (C) is a multibody analysis code.\n');
fprintf(fd, '# http://www.mbdyn.org\n');
fprintf(fd, '# \n');
fprintf(fd, '# Copyright (C) 1996-2021\n');
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
fprintf(fd, '\n');

%m = 11.64;
%JXX = 13.17e-3;
%J11 = 2*JXX;
mD = 70.573; % kg
rD = 0.24; % m
tD = 0.05; % m
JD11 = 1e-0*mD*rD^2/2;
JDXX = 1e-0*mD*(rD^2/4 + tD^2/12);
%JDXX = 1.0163e-3;
%JD11 = 2*JDXX;
D = 0.05; %

for i = 0:2*N_BEAMS,
	DL = L/(2*N_BEAMS);
	DX = 0.;
	if (i == 0),
		DL = DL/2;
		DX = DL/4;
	elseif (i == 2*N_BEAMS),
		DL = DL/2;
		DX = -DL/4;
	end
	DM = DL*m;
	DJ11 = DL*J11;
	DJXX = DM*DL^2/12 + DL*JXX;

	fprintf(fd, 'body: %d, %d', i, i);
	if (i == N_BEAMS),
		fprintf(fd, ',\n');
		fprintf(fd, '\tcondense, 2');
	end
	fprintf(fd, ',\n');
	fprintf(fd, '\t%+20.10e,\n', DM);
	fprintf(fd, '\treference, node, %+20.10e, %+20.10e, %+20.10e,\n', DX, 0., 0.);
	fprintf(fd, '\tdiag, %+20.10e, %+20.10e, %+20.10e', DJ11, DJXX, DJXX);
	if (i == N_BEAMS),
		fprintf(fd, ',\n');
		fprintf(fd, '\t%+20.10e,\n', mD);
		fprintf(fd, '\treference, node, %+20.10e, %+20.10e, %+20.10e,\n', 0., 0., D);
		fprintf(fd, '\tdiag, %+20.10e, %+20.10e, %+20.10e', JD11, JDXX, JDXX);
	end
	fprintf(fd, ';\n');

	fprintf(fd, '\n');
end

for i = 1:N_BEAMS,
	fprintf(fd, 'beam3: %d,\n', i);
	for j = 1:3,
		fprintf(fd, '\t%d,\n', 2*i + j - 3);
		fprintf(fd, '\t\tposition, reference, node, null,\n');
		fprintf(fd, '\t\torientation, reference, node, eye,\n');
	end
	fprintf(fd, '\tfrom nodes,\n');
	fprintf(fd, '\treference, 6,\n');
	fprintf(fd, '\tsame,\n');
	fprintf(fd, '\tsame;\n');

	fprintf(fd, '\n');
end

fclose(fd);

