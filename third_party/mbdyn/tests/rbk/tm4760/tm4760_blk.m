% $Header: /var/cvs/mbdyn/mbdyn/mbdyn-1.0/tests/rbk/tm4760/tm4760_blk.m,v 1.5 2008/12/31 14:02:09 masarati Exp $
%
% NASA Technical Memorandum 4760
% ARL Technical Report 1389
%
% Rotating Shake Test and Modal Analysis of a Model Helicopter Rotor Blade
% W. Keats Wilkie, Paul H. Mirick, and Chester W. Langston
% Vehicle Technology Center
% U.S. Army Research Laboratory
% Langley Research Center Hampton, Virginia
%
% Data from 'tm4760.blk'

% beamtype:
%	2: use beam2 (shear stiffness modified to eliminate shear locking)
%	3: use beam3 (adding intermediate static node)
beamtype = 2;

% GRID
GRID = [...
% 1-GID           2-X  
%                 in
204               3.00
200               3.00
201               6.87
202               8.87
203              10.625
%
1000             12.500
1001             13.000
1002             14.000
1003             15.000
1004             15.375
1005             15.600
1006             16.000
1007             17.000
1008             17.850
1009             18.000
1010             19.000
1011             20.000
1012             21.000
1013             22.000
1014             23.000
1015             23.750
1016             24.000
1017             25.000
1018             26.000
1019             27.000
1020             28.000
1021             28.250
1022             29.000
1023             29.150
1024             30.000
1025             31.000
1026             32.000
1027             33.000
1028             34.000
1029             35.000
1030             36.000
1031             37.000
1032             38.000
1033             39.000
1034             40.000
1035             41.000
1036             41.250
1037             42.000
1038             43.000
1039             44.000
1040             45.000
1041             46.000
1042             47.000
1043             48.000
1044             49.000
1045             50.000
1046             51.000
1047             52.750
1048             53.000
1049             54.000
1050             54.250
1051             55.000
];

% CBEAM
CBEAM = [...
% 1-CID 2-PID   3-G1    4-G2
200     200     200     201
201     201     201     202
202     202     202     203
203     203     203     1000
%
1000    1000    1000    1001
1001    1001    1001    1002
1002    1001    1002    1003
1003    1001    1003    1004
1004    1004    1004    1005
1005    1004    1005    1006
1006    1004    1006    1007
1007    1004    1007    1008
1008    1008    1008    1009
1009    1008    1009    1010
1010    1008    1010    1011
1011    1008    1011    1012
1012    1008    1012    1013
1013    1008    1013    1014
1014    1008    1014    1015
1015    1015    1015    1016
1016    1015    1016    1017
1017    1015    1017    1018
1018    1015    1018    1019
1019    1015    1019    1020
1020    1015    1020    1021
1021    1021    1021    1022
1022    1021    1022    1023
1023    1021    1023    1024
1024    1021    1024    1025
1025    1021    1025    1026
1026    1021    1026    1027
1027    1021    1027    1028
1028    1021    1028    1029
1029    1021    1029    1030
1030    1021    1030    1031
1031    1021    1031    1032
1032    1021    1032    1033
1033    1021    1033    1034
1034    1021    1034    1035
1035    1021    1035    1036
1036    1021    1036    1037
1037    1021    1037    1038
1038    1021    1038    1039
1039    1021    1039    1040
1040    1021    1040    1041
1041    1021    1041    1042
1042    1021    1042    1043
1043    1021    1043    1044
1044    1021    1044    1045
1045    1021    1045    1046
1046    1046    1046    1047
1047    1047    1047    1048
1048    1048    1048    1049
1049    1049    1049    1050
1050    1050    1050    1051
];

% PBEAM
PBEAM = [...
% 1-PID 2-A     3-I1    4-I2    5-IP    6-M       7-JP      8-YCM
%       in^2    in^4    in^4    in^4    lb/in     lb-in^2/in in
200     5.500   0.5000  0.5000  0.2632  0.42511   2.200e-1  0.00
201     0.371   0.1500  0.0500  0.1316  0.19380   4.826e-2  0.00
202     0.371   0.0250  0.0040  0.0105  0.04086   2.505e-2  0.00
203     0.371   0.0250  0.0040  0.0105  0.15113   5.671e-2  0.00
1000    0.377   0.03550 0.00394 0.03940 0.14048   5.559e-2  0.00
1001    0.386   0.02520 0.00249 0.00976 0.3134e-1 2.779e-2  0.00
1004    0.339   0.02520 0.00249 0.00976 0.4376e-1 2.961e-2  0.00
1008    0.278   0.03040 0.00231 0.00674 0.4123e-1 2.888e-2  0.00
1015    0.249   0.02636 0.00181 0.00565 0.400e-1  2.810e-2  0.00
1021    0.224   0.02447 0.00151 0.00485 0.3903e-1 2.756e-2  0.00
1046    0.247   0.02447 0.00151 0.00485 0.3903e-1 2.756e-2  0.00
1047    0.279   0.02448 0.00160 0.00502 0.4125e-1 2.826e-2  0.00
1048    0.305   0.05000 0.00500 0.01435 0.7867e-1 3.984e-2 -0.106
1049    0.099   0.04000 0.00400 0.01148 0.607e-1  3.563e-2 -0.17
1050    0.05    0.00500 0.00050 0.01435 0.1101e-1 4.401e-3 -0.940
];

[ngrid, dmy] = size(GRID);
[ncbeam, dmy] = size(CBEAM);
[npbeam, dmy] = size(PBEAM);

%%%%%%%%%%%
%% setup %%
%%%%%%%%%%%
fd = fopen('tm4760_blk.set', 'w');

fprintf(fd, '# $Header: /var/cvs/mbdyn/mbdyn/mbdyn-1.0/tests/rbk/tm4760/tm4760_blk.m,v 1.5 2008/12/31 14:02:09 masarati Exp $\n');
fprintf(fd, '#\n');
fprintf(fd, '# Rotor blade model for MBDyn > 1.3.4\n');
fprintf(fd, '# Author: Pierangelo Masarati <masarati@aero.polimi.it>\n');
fprintf(fd, '# Data from:\n');
fprintf(fd, '#\n');
fprintf(fd, '# NASA Technical Memorandum 4760\n');
fprintf(fd, '# ARL Technical Report 1389\n');
fprintf(fd, '#\n');
fprintf(fd, '# Rotating Shake Test and Modal Analysis of a Model Helicopter Rotor Blade\n');
fprintf(fd, '# W. Keats Wilkie, Paul H. Mirick, and Chester W. Langston\n');
fprintf(fd, '# Vehicle Technology Center\n');
fprintf(fd, '# U.S. Army Research Laboratory\n');
fprintf(fd, '# Langley Research Center Hampton, Virginia\n');
fprintf(fd, '#\n');
fprintf(fd, '# file generated by ''tm4760_blk.m'' - modifications could be overwritten\n');
fprintf(fd, '# rotor blade setup for MBDyn > 1.3.4\n');

if (beamtype == 2),
	NSTRUCTNODE = ngrid;
else
	NSTRUCTNODE = ngrid + ncbeam;
end
NBODY = ngrid - 1;
NBEAM = ncbeam;

fprintf(fd, '\n');
fprintf(fd, 'set: integer NSTRUCTNODE = %d;\n', NSTRUCTNODE);
fprintf(fd, 'set: integer NBODY = %d;\n', NBODY);
fprintf(fd, 'set: integer NBEAM = %d;\n', NBEAM);

fclose(fd);

%%%%%%%%%%%
%% nodes %%
%%%%%%%%%%%
fd = fopen('tm4760_blk.nod', 'w');

fprintf(fd, '# $Header: /var/cvs/mbdyn/mbdyn/mbdyn-1.0/tests/rbk/tm4760/tm4760_blk.m,v 1.5 2008/12/31 14:02:09 masarati Exp $\n');
fprintf(fd, '#\n');
fprintf(fd, '# Rotor blade model for MBDyn > 1.3.4\n');
fprintf(fd, '# Author: Pierangelo Masarati <masarati@aero.polimi.it>\n');
fprintf(fd, '# Data from:\n');
fprintf(fd, '#\n');
fprintf(fd, '# NASA Technical Memorandum 4760\n');
fprintf(fd, '# ARL Technical Report 1389\n');
fprintf(fd, '#\n');
fprintf(fd, '# Rotating Shake Test and Modal Analysis of a Model Helicopter Rotor Blade\n');
fprintf(fd, '# W. Keats Wilkie, Paul H. Mirick, and Chester W. Langston\n');
fprintf(fd, '# Vehicle Technology Center\n');
fprintf(fd, '# U.S. Army Research Laboratory\n');
fprintf(fd, '# Langley Research Center Hampton, Virginia\n');
fprintf(fd, '#\n');
fprintf(fd, '# file generated by ''tm4760_blk.m'' - modifications could be overwritten\n');
fprintf(fd, '# rotor blade nodes for MBDyn > 1.3.4\n');

%%%%%%%%%%%%%%%%%%%%%%
%% structural nodes %%
%%%%%%%%%%%%%%%%%%%%%%

fprintf(fd, '\n');
fprintf(fd, '# structural nodes\n');

% first node, clamped
fprintf(fd, '\n');
fprintf(fd, 'set: integer GRID_%d = %d;\n', GRID(1, 1), GRID(1, 1));
fprintf(fd, 'structural: GRID_%d, static,\n', GRID(1, 1));
fprintf(fd, '\tposition, reference, global,\n');
fprintf(fd, '\t\t%+17.10e, %+17.10e, %+17.10e,\n', GRID(1, 2), 0., 0.);
fprintf(fd, '\torientation, reference, global, eye,\n');
fprintf(fd, '\tvelocity, reference, global, null,\n');
fprintf(fd, '\tangular velocity, reference, global, null;\n');

for n = 2:ngrid,
	if (beamtype == 3 && n > 2),
		fprintf(fd, '\n');
		fprintf(fd, 'set: integer GRID_%d_MID = %d;\n', GRID(n, 1), 100000 + GRID(n, 1));
		fprintf(fd, 'structural: GRID_%d_MID, static,\n', GRID(n, 1));
		fprintf(fd, '\tposition, reference, global,\n');
		fprintf(fd, '\t\t%+17.10e, %+17.10e, %+17.10e,\n', (GRID(n - 1, 2) + GRID(n, 2))/2., 0., 0.);
		fprintf(fd, '\torientation, reference, global, eye,\n');
		fprintf(fd, '\tvelocity, reference, global, null,\n');
		fprintf(fd, '\tangular velocity, reference, global, null;\n');
	end

	fprintf(fd, '\n');
	fprintf(fd, 'set: integer GRID_%d = %d;\n', GRID(n, 1), GRID(n, 1));
	fprintf(fd, 'structural: GRID_%d, dynamic,\n', GRID(n, 1));
	fprintf(fd, '\tposition, reference, global,\n');
	fprintf(fd, '\t\t%+17.10e, %+17.10e, %+17.10e,\n', GRID(n, 2), 0., 0.);
	fprintf(fd, '\torientation, reference, global, eye,\n');
	fprintf(fd, '\tvelocity, reference, global, null,\n');
	fprintf(fd, '\tangular velocity, reference, global, null;\n');
end

fclose(fd);

%%%%%%%%%%%%%%
%% elements %%
%%%%%%%%%%%%%%
fd = fopen('tm4760_blk.elm', 'w');

fprintf(fd, '# $Header: /var/cvs/mbdyn/mbdyn/mbdyn-1.0/tests/rbk/tm4760/tm4760_blk.m,v 1.5 2008/12/31 14:02:09 masarati Exp $\n');
fprintf(fd, '#\n');
fprintf(fd, '# Rotor blade model for MBDyn > 1.3.4\n');
fprintf(fd, '# Author: Pierangelo Masarati <masarati@aero.polimi.it>\n');
fprintf(fd, '# Data from:\n');
fprintf(fd, '#\n');
fprintf(fd, '# NASA Technical Memorandum 4760\n');
fprintf(fd, '# ARL Technical Report 1389\n');
fprintf(fd, '#\n');
fprintf(fd, '# Rotating Shake Test and Modal Analysis of a Model Helicopter Rotor Blade\n');
fprintf(fd, '# W. Keats Wilkie, Paul H. Mirick, and Chester W. Langston\n');
fprintf(fd, '# Vehicle Technology Center\n');
fprintf(fd, '# U.S. Army Research Laboratory\n');
fprintf(fd, '# Langley Research Center Hampton, Virginia\n');
fprintf(fd, '#\n');
fprintf(fd, '# file generated by ''tm4760_blk.m'' - modifications could be overwritten\n');
fprintf(fd, '# rotor blade elements for MBDyn > 1.3.4\n');

%%%%%%%%%%%%%%%%%%
%% rigid bodies %%
%%%%%%%%%%%%%%%%%%

fprintf(fd, '\n');
fprintf(fd, '# rigid bodies\n');

% first GRID is clamped and serves as pin for the hinge
for b = 2:ngrid,
	if (b == 2 || b == ngrid),
		condense = 1;
	else
		condense = 2;
	end

	fprintf(fd, '\n');
	fprintf(fd, 'body: GRID_%d, GRID_%d', GRID(b, 1), GRID(b, 1));
	if (condense > 1),
		fprintf(fd, ',\n\tcondense, %d', condense);
	end

	% CBEAM node b is g2 of
	cbeam_idx = find(CBEAM(:, 4) == GRID(b, 1));
	if (isempty(cbeam_idx) == 0),
		g1_idx = find(GRID(:, 1) == CBEAM(cbeam_idx, 3));
		if (isempty(g1_idx)),
			error(sprintf('body(%d)=%d: unable to find G1=%d of CBEAM %d', ...
				b, GRID(b, 1), CBEAM(cbeam_idx, 3), CBEAM(cbeam_idx, 1)));
		end

		pbeam_idx = find(PBEAM(:, 1) == CBEAM(cbeam_idx, 2));
		if (isempty(pbeam_idx)),
			error(sprintf('body(%d)=%d: unable to find PBEAM=%d of CBEAM %d', ...
				b, GRID(b, 1), CBEAM(cbeam_idx, 2), CBEAM(cbeam_idx, 1)));
		end

		dl = GRID(b, 2) - GRID(g1_idx, 2);
		dm = PBEAM(pbeam_idx, 6);
		dJp = PBEAM(pbeam_idx, 7);
		ycm = PBEAM(pbeam_idx, 8);

		fprintf(fd, ',\n\t# CBEAM=%d, PBEAM=%d, G1=%d, G2=%d, L=%.4f in\n', ...
			CBEAM(cbeam_idx, 1), CBEAM(cbeam_idx, 2), CBEAM(cbeam_idx, 3), CBEAM(cbeam_idx, 4), dl/2.);
		fprintf(fd, '\t%+17.10e,\n', dm*dl/2.);
		fprintf(fd, '\treference, node, %+17.10e, %+17.10e, %+17.10e,\n', -dl/4., ycm, 0.);
		fprintf(fd, '\tdiag, %+17.10e, %+17.10e, %+17.10e', dJp*dl/2., dm/12*(dl/2.)^3, dm/12*(dl/2.)^3);
	end

	% CBEAM node b is g1 of
	cbeam_idx = find(CBEAM(:, 3) == GRID(b, 1));
	if (isempty(cbeam_idx) == 0),
		g2_idx = find(GRID(:, 1) == CBEAM(cbeam_idx, 4));
		if (isempty(g2_idx)),
			error(sprintf('body(%d)=%d: unable to find G2=%d of CBEAM %d', ...
				b, GRID(b, 1), CBEAM(cbeam_idx, 4), CBEAM(cbeam_idx, 1)));
		end

		pbeam_idx = find(PBEAM(:, 1) == CBEAM(cbeam_idx, 2));
		if (isempty(pbeam_idx)),
			error(sprintf('body(%d)=%d: unable to find PBEAM=%d of CBEAM %d', ...
				b, GRID(b, 1), CBEAM(cbeam_idx, 2), CBEAM(cbeam_idx, 1)));
		end

		dl = GRID(g2_idx, 2) - GRID(b, 2);
		dm = PBEAM(pbeam_idx, 6);
		dJp = PBEAM(pbeam_idx, 7);
		ycm = PBEAM(pbeam_idx, 8);

		fprintf(fd, ',\n\t# CBEAM=%d, PBEAM=%d, G1=%d, G2=%d, L=%.4f in\n', ...
			CBEAM(cbeam_idx, 1), CBEAM(cbeam_idx, 2), CBEAM(cbeam_idx, 3), CBEAM(cbeam_idx, 4), dl/2.);
		fprintf(fd, '\t%+17.10e,\n', dm*dl/2.);
		fprintf(fd, '\treference, node, %+17.10e, %+17.10e, %+17.10e,\n', dl/4., ycm, 0.);
		fprintf(fd, '\tdiag, %+17.10e, %+17.10e, %+17.10e', dJp*dl/2., dm/12*(dl/2.)^3, dm/12*(dl/2.)^3);
	end

	fprintf(fd, ';\n');
end

fprintf(fd, '\n');
fprintf(fd, 'inertia: 0, body, all;\n');

%%%%%%%%%%%
%% beams %%
%%%%%%%%%%%

fprintf(fd, '\n');
fprintf(fd, '# beam constitutive laws\n');

fprintf(fd, '\n');
fprintf(fd, 'set: real E = %+7.10e*(9.81/.3048*12);\t# lb/in^2\n', 1e+7);
fprintf(fd, 'set: real NU = %+7.10e;\n', 0.3);
fprintf(fd, 'set: real G = E/(2*(1 + NU));\n');
fprintf(fd, 'set: real F_DAMP = %+17.10e;\n', 0.);

if (beamtype == 3),
	for p = 1:npbeam,
		fprintf(fd, '\n');
		fprintf(fd, 'set: integer PBEAM_%d = %d;\n', PBEAM(p, 1), PBEAM(p, 1));
		fprintf(fd, 'constitutive law: PBEAM_%d, 6,\n', PBEAM(p, 1));
		fprintf(fd, '\tlinear viscoelastic generic, diag,\n');
		fprintf(fd, '\t\t%16.10e*E,\n', PBEAM(p, 2));
		fprintf(fd, '\t\t%16.10e*G,\n', PBEAM(p, 2));
		fprintf(fd, '\t\t%16.10e*G,\n', PBEAM(p, 2));
		fprintf(fd, '\t\t%16.10e*G,\n', PBEAM(p, 5));
		fprintf(fd, '\t\t%16.10e*E,\n', PBEAM(p, 4));
		fprintf(fd, '\t\t%16.10e*E,\n', PBEAM(p, 3));
		fprintf(fd, '\tproportional, F_DAMP;\n');
	end

else
	for b = 1:ncbeam,
		pbeam_idx = find(PBEAM(:, 1) == CBEAM(b, 2));
		g1_idx = find(GRID(:, 1) == CBEAM(b, 3));
		g2_idx = find(GRID(:, 1) == CBEAM(b, 4));

		L = GRID(g2_idx, 2) - GRID(g1_idx, 2);
		A = PBEAM(pbeam_idx, 2);
		I1 = PBEAM(pbeam_idx, 3);
		I2 = PBEAM(pbeam_idx, 4);
		IP = PBEAM(pbeam_idx, 5);

		fprintf(fd, '\n');
		fprintf(fd, 'set: integer PBEAM_%d = %d;\n', CBEAM(b, 1), CBEAM(b, 1));
		fprintf(fd, 'constitutive law: PBEAM_%d, 6,\n', CBEAM(b, 1));
		fprintf(fd, '\tlinear viscoelastic generic, diag,\n');
		fprintf(fd, '\t\t%16.10e*E,\n', A);
		fprintf(fd, '\t\t1./(1./(%16.10e*G) + %16.10e^2/(12.*%16.10e*E)),\n', A, L, I1);
		fprintf(fd, '\t\t1./(1./(%16.10e*G) + %16.10e^2/(12.*%16.10e*E)),\n', A, L, I2);
		fprintf(fd, '\t\t%16.10e*G,\n', IP);
		fprintf(fd, '\t\t%16.10e*E,\n', I2);
		fprintf(fd, '\t\t%16.10e*E,\n', I1);
		fprintf(fd, '\tproportional, F_DAMP;\n');
	end
end

fprintf(fd, '\n');
fprintf(fd, '# beams\n');

if (beamtype == 3),
	for b = 1:ncbeam,
		fprintf(fd, '\n');
		fprintf(fd, 'beam3: %d,\n', CBEAM(b, 1));
		fprintf(fd, '\tGRID_%d, null,\n', CBEAM(b, 3));
		fprintf(fd, '\tGRID_%d_MID, null,\n', CBEAM(b, 4));
		fprintf(fd, '\tGRID_%d, null,\n', CBEAM(b, 4));
		fprintf(fd, '\treference, global, eye,\n');
		fprintf(fd, '\treference, PBEAM_%d,\n', CBEAM(b, 2));
		fprintf(fd, '\tsame,\n');
		fprintf(fd, '\tsame;\n');
	end

else
	for b = 1:ncbeam,
		fprintf(fd, '\n');
		fprintf(fd, 'beam2: %d,\n', CBEAM(b, 1));
		fprintf(fd, '\tGRID_%d, null,\n', CBEAM(b, 3));
		fprintf(fd, '\tGRID_%d, null,\n', CBEAM(b, 4));
		fprintf(fd, '\treference, global, eye,\n');
		fprintf(fd, '\treference, PBEAM_%d;\n', CBEAM(b, 1));
	end
end

fclose(fd);
