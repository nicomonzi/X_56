% $Header: /var/cvs/mbdyn/mbdyn/mbdyn-1.0/tests/forces/modalext/socket/modalmapping.m,v 1.1 2009/12/14 21:43:12 masarati Exp $

CHORD = .1;

NNODES = 2;
XI = [0:NNODES]'/NNODES;

XB = [0; 1; 2]/2;
XC = [0; 1; 4]/4;
XT = [1; 1; 1];

NMODESB = 1;
NMODESC = 1;
NMODEST = 1;
H = zeros(6*(NNODES + 1), NMODESB + NMODESC + NMODEST);
M = zeros(6*(NNODES + 1), 6*(NNODES + 1));
for i = 1:NNODES + 1,
	% chordwise caused by bending
	H(6*(i - 1) + 2, NMODESB + [1:NMODESC]) = XC(i, 1:NMODESC);
	% beamwise caused by bending
	H(6*(i - 1) + 3, 0 + [1:NMODESB]) = XB(i, 1:NMODESB);
	% beamwise caused by torsion
	H(6*(i - 1) + 3, NMODESB + NMODESC + [1:NMODEST]) = .5*CHORD*XT(i, 1:NMODEST);

	% chordwise caused by bending
	H(6*(i - 1) + 5, NMODESB + [1:NMODESC]) = XC(i, 1:NMODESC);
	% beamwise caused by bending
	H(6*(i - 1) + 6, 0 + [1:NMODESB]) = XB(i, 1:NMODESB);
	% beamwise caused by torsion
	H(6*(i - 1) + 6, NMODESB + NMODESC + [1:NMODEST]) = -.5*CHORD*XT(i, 1:NMODEST);

	% leading edge
	M(6*(i - 1) + [1:3], 6*(i - 1) + [1:3]) = eye(3);
	l = [0.; .5*CHORD; 0.];
	M(6*(i - 1) + [1:3], 6*(i - 1) + [4:6]) = -[0, -l(3), l(2); l(3), 0, -l(1); -l(2), l(1), 0];

	% trailing edge
	M(6*(i - 1) + [4:6], 6*(i - 1) + [1:3]) = eye(3);
	t = [0.; -.5*CHORD; 0.];
	M(6*(i - 1) + [4:6], 6*(i - 1) + [4:6]) = -[0, -t(3), t(2); t(3), 0, -t(1); -t(2), t(1), 0];
end

Hp = H\M;

Hpx = Hp;
IDX = find(abs(Hp) < 1e-13);
Hpx(IDX) = 0;

Hps = sparse(Hpx);

save 'Hp_full.dat' Hpx
save 'Hp_sparse.dat' Hps

