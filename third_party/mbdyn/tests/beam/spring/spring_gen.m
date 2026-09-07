% $Header: /var/cvs/mbdyn/mbdyn/mbdyn-1.0/tests/beam/spring/spring.m,v 1.2 2007/02/26 00:34:35 masarati Exp $

% helix
nbeams = 50;
L = 0.5;	% m
p = 0.1;
R0 = .05;
Rmax = .08;

EA = 1e+8;
GA = 1e+7;
GJ = 1e+4;
EJ = 1e+4;
blade_damp = 1.e-4;

t = [-nbeams:nbeams]/nbeams;
R = Rmax - (Rmax - R0)*t.^2;
x = cos(t*pi*L/p).*R;
y = sin(t*pi*L/p).*R;
z = [t]*L/2;

xp = -pi*L/p*sin(t*pi*L/p).*R - 2*(Rmax - R0)*cos(t*pi*L/p).*t;
yp = pi*L/p*cos(t*pi*L/p).*R - 2*(Rmax - R0)*sin(t*pi*L/p).*t;
zp = ones(size(t))*L/2;

fd = fopen('spring.nod', 'w');

fprintf(fd, '# ground\n');
fprintf(fd, 'structural: 0, static,\n');
fprintf(fd, '        null,\n');
fprintf(fd, '        eye,\n');
fprintf(fd, '        null,\n');
fprintf(fd, '        null;\n');
fprintf(fd, '\n');

fprintf(fd, '# top\n');
fprintf(fd, 'structural: spring_top, dynamic,\n');
fprintf(fd, '        reference, spring_top, null,\n');
fprintf(fd, '        reference, spring_top, eye,\n');
fprintf(fd, '        reference, spring_top, null,\n');
fprintf(fd, '        reference, spring_top, null;\n');
fprintf(fd, '\n');

for i = 0:2*nbeams,
	fprintf(fd, '# node %d\n', i);
	fprintf(fd, 'structural: spring + %d, dynamic,\n', i);
	fprintf(fd, '        reference, spring, %+e,%+e,%+e,\n', ...
		x(i + 1), y(i + 1), z(i + 1));
	fprintf(fd, '        reference, spring,\n');
	fprintf(fd, '                1, %+e,%+e,%+e,\n', ...
		xp(i + 1), yp(i + 1), zp(i + 1));
	fprintf(fd, '                2, %+e,%+e,%+e,\n', ...
		-yp(i + 1), xp(i + 1), 0.);
	fprintf(fd, '        reference, spring, null,\n');
	fprintf(fd, '        reference, spring, null;\n');
	fprintf(fd, '\n');
end

fclose(fd);

fopen('spring.elm', 'w');

fprintf(fd, '# clamp for ground\n');
fprintf(fd, 'joint: spring + 1000, clamp, 0, node, node;\n');
fprintf(fd, '\n');

fprintf(fd, '# inline for top\n');
fprintf(fd, 'joint: spring + 1001, inline,\n');
fprintf(fd, '        0,\n');
fprintf(fd, '                reference, spring_top, null,\n');
fprintf(fd, '                reference, spring_top, eye,\n');
fprintf(fd, '        spring_top;\n');
fprintf(fd, '\n');

fprintf(fd, '# prismatic for top\n');
fprintf(fd, 'joint: spring + 1002, prismatic,\n');
fprintf(fd, '        0,\n');
fprintf(fd, '        spring_top;\n');
fprintf(fd, '\n');

fprintf(fd, '# body for top\n');
fprintf(fd, 'body: spring_top, spring_top,\n');
fprintf(fd, '        100.,\n');
fprintf(fd, '        null,\n');
fprintf(fd, '        diag, 1.,1.,1.;\n');
fprintf(fd, '\n');

fprintf(fd, '# clamp at bottom\n');
fprintf(fd, 'joint: spring + 1003, clamp, spring + 0, node, node;\n');
fprintf(fd, '\n');

fprintf(fd, '# inplane at top\n');
fprintf(fd, 'joint: spring + 1004, inplane,\n');
fprintf(fd, '        spring_top,\n');
fprintf(fd, '                reference, spring_top, null,\n');
fprintf(fd, '                reference, spring_top, 0.,0.,1.,\n');
fprintf(fd, '        spring + %d;\n', 2*nbeams);

for i = 0:2*nbeams,
	fprintf(fd, '# body %d\n', i);
	fprintf(fd, 'body: spring + %d, spring + %d,\n', i, i);
	fprintf(fd, '        %+e,\n', 1.e-3);
	fprintf(fd, '        null,\n');
	fprintf(fd, '        diag, %+e,%+e,%+e;\n', 1.e-4, 1.e-3, 1.e-3);
	fprintf(fd, '\n');
end

tb(1:2:2*nbeams) = [-nbeams + (1 - 1/sqrt(3)):2:nbeams]/nbeams;
tb(2:2:2*nbeams) = [-nbeams + (1 + 1/sqrt(3)):2:nbeams]/nbeams;

Rb = Rmax - (Rmax - R0)*tb.^2;
xpb = -pi*L/p*sin(tb*pi*L/p).*Rb - 2*(Rmax - R0)*cos(tb*pi*L/p).*tb;
ypb = pi*L/p*cos(tb*pi*L/p).*Rb - 2*(Rmax - R0)*sin(tb*pi*L/p).*tb;
zpb = ones(size(tb))*L/2;

for i = 1:nbeams,
	fprintf(fd, '# beam %d\n', i);
	fprintf(fd, 'beam3: spring + %d,\n', i);
	fprintf(fd, '        spring + %d, null,\n', 2*i - 2);
	fprintf(fd, '        spring + %d, null,\n', 2*i - 1);
	fprintf(fd, '        spring + %d, null,\n', 2*i - 0);
	fprintf(fd, '        reference, spring,\n');
	fprintf(fd, '                1, %+e,%+e,%+e,\n', ...
		xpb(2*i - 1), ypb(2*i - 1), zpb(2*i - 1));
	fprintf(fd, '                2, %+e,%+e,%+e,\n', ...
		xpb(2*i), ypb(2*i), zpb(2*i));
	fprintf(fd, '        linear viscoelastic generic, diag,\n');
	fprintf(fd, '                %+e,%+e,%+e,\n', ...
		EA, GA, GA);
	fprintf(fd, '                %+e,%+e,%+e,\n', ...
		GJ, EJ, EJ);
	fprintf(fd, '                proportional, %+e,\n', blade_damp);
	fprintf(fd, '        reference, spring,\n');
	fprintf(fd, '                1, %+e,%+e,%+e,\n', ...
		xpb(2*i - 1), ypb(2*i - 1), zpb(2*i - 1));
	fprintf(fd, '                2, %+e,%+e,%+e,\n', ...
		xpb(2*i), ypb(2*i), zpb(2*i));
	fprintf(fd, '        same;\n');
	fprintf(fd, '\n');
end

fclose(fd);
