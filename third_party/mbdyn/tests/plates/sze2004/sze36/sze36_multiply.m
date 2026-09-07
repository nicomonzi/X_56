% file zz.mov contains .mov at step 0 and at the desired step for output
mov = load('zz.mov');

%
%    z            +------+------+
%    ^            |      |      |
%    |            | 1, 2 | 3, 4 |
%    |            |      |      |
%    x ---> y     +------+------+
%                 |      |      |
%                 | 5, 6 | 7, 8 |
%                 |      |      |
%                 +------+------+
%

signs = [1,1,1; -1,1,1; 1,-1,1; -1,-1,1; 1,1,-1; -1,1,-1; 1,-1,-1; -1,-1,-1];
mov = mov(:, 1:4);
offsets = [0,0,0; 0,0,0; 0,200,0; 0,200,0; 0,0,0; 0,0,0; 0,200,0; 0,200,0];

nn = size(mov, 1)/2;
newmov = zeros(8*nn*2, 4);

for t = 1:2,
	newmov((t-1)*8*nn + 1:((t-1)*8 + 1)*nn, :) = mov((t-1)*nn+1:t*nn, :);
	for b = 2:8,
		% label
		newmov(((t-1)*8 + b - 1)*nn + 1:((t-1)*8 + b)*nn, 1) = ...
			(b-1)*1000000 + mov((t-1)*nn+1:t*nn, 1);
		for c = 1:3,
			newmov(((t-1)*8 + b - 1)*nn + 1:((t-1)*8 + b)*nn, 1 + c) = ...
				offsets(b, c) + signs(b, c)*mov((t-1)*nn+1:t*nn, 1 + c);
		end
	end
end

fid = fopen('z2.mov', 'w');
n = size(newmov, 1);
for i = 1:n,
	fprintf(fid, '%d %e %e %e 0 0 0 0 0 0 0 0 0\n', ...
		newmov(i, 1), newmov(i, 2), newmov(i, 3), newmov(i, 4));
end
fclose(fid);

fid = fopen('z2.log', 'w');

% file zz.strnode contains "grep '^structural node:' | sed -e 's/^structural node: //' -e 's/euler123 //'" from .log
nod = load('zz.strnode');
n = size(nod, 1);
for b = 1:8,
	for i = 1:n,
		fprintf(fid, 'structural node: %d %e %e %e euler123 0 0 0\n', ...
			newmov((b-1)*n + i, 1), ...
			newmov((b-1)*n + i, 2), ...
			newmov((b-1)*n + i, 3), ...
			newmov((b-1)*n + i, 4));
	end
end

% file zz.shell contains "grep '^shell4:' | sed -e 's/^shell4: //'" from .log
she = load('zz.shell');
n = size(she, 1);
for b = 1:8,
	if (signs(b, 1)*signs(b, 2)*signs(b, 3) < 0),
		n1 = 5;
		n2 = 4;
		n3 = 3;
		n4 = 2;
	else
		n1 = 2;
		n2 = 3;
		n3 = 4;
		n4 = 5;
	end
	for i = 1:n,
		fprintf(fid, 'shell4: %d %d %d %d %d\n', ...
			(b-1)*1000000 + she(i, 1), ...
			(b-1)*1000000 + she(i, n1), ...
			(b-1)*1000000 + she(i, n2), ...
			(b-1)*1000000 + she(i, n3), ...
			(b-1)*1000000 + she(i, n4));

	end
end

fclose(fid);
