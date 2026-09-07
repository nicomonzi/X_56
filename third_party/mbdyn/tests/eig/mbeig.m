% $Header: /var/cvs/mbdyn/mbdyn/mbdyn-1.0/tests/eig/mbeig.m,v 1.13 2011/08/26 01:08:08 masarati Exp $
%
% If "output matrices" is set, eigenvalues and eigenvectors are computed
% by octave/matlab, and the reduced model is generated.
%
% If "ouput eigenvectors" is set, eigenvalues and eigenvectors computed
% by MBDyn are used; only if "output matrices" is set, the reduced model
% is generated.
%
% NOTE: you may set "upper frequency limit", but you should not
% set "lower frequency limit"

% set frequency range of interest
if exist('fmin') ~= 1,
	fmin = 1.e-3*2*pi;
end
if exist('fmax') ~= 1,
	fmax = 1.e+2*2*pi;
end

if exist('VR') == 1 & exist('alpha') == 1,
	if (isempty(find(abs(alpha(:, 3)) < 1e-12)) ~= 1),
		error('Houston, we''ve had a problem: (alpha_r + j*alpha_i)/beta: some beta are too small');
	end
	Y = VR;
	L = (alpha(:, 1) + sqrt(-1)*alpha(:, 2))./alpha(:, 3);

	if (exist('VL') == 1),
		Yc = VL;
	end

else
	% "right" and "left" eigenproblem
	if exist('OCTAVE_HOME'),
		[dmy, dmy, dmy, dmy, Y, Yc, L] = qz(Aminus, Aplus);
		[n, dmy] = size(Aminus);
		nL = length(L);
		if (nL < n),
			error(sprintf('Houston, we''ve had a problem: too few valid eigenvalues (%d instead of %d)...', nL, n));
		end
	else
		[AA, BB, dmy, dmy, Y, Yc] = qz(Aminus, Aplus);
		% NOTE: by default, matlab performs the complex QZ transform.
		% As a consequence, the eigenvalues are diag(AA)./diag(BB).
		% This behavior is indicated as legacy; it might change
		% in the future.
		L = diag(AA)./diag(BB);
	end
end

% reduce "right" problem
l = (L - 1)./(L + 1)/dCoef;
if exist('vv') == 1,
	v = vv;
else
	v = find((abs(imag(l)) >= fmin) & (abs(imag(l)) <= fmax));
end

m = length(v);

% combine complex modes; leave real alone
T = zeros(m);
i = 1;
while (i <= m),
	if (norm(imag(Y(:, v(i)))) == 0.),
		T(i, i) = 1;
		i = i + 1;
	else
		T(i:i + 1, i:i + 1) = [.5 -.5*sqrt(-1); .5 .5*sqrt(-1)];
		i = i + 2;
	end
end
y = Y(:, v)*T;

% need left eigenvectors for model reduction
if (exist('Yc') ~= 1),
	return;
end

% matrices Aplus, Aminus and dCoef from MBDyn's run must be available
if exist('Aplus') ~= 1 | exist('Aminus') ~= 1 | exist('dCoef') ~= 1,
	error('need (at least) Aplus, Aminus and dCoef; execute ''<outfile>.m'' first');
end

Tc = zeros(m);
i = 1;
while (i <= m),
	if (norm(imag(Yc(:, v(i)))) == 0.),
		Tc(i, i) = 1;
		i = i + 1;
	else
		Tc(i:i + 1, i:i + 1) = [.5 -.5*sqrt(-1); .5 .5*sqrt(-1)];
		i = i + 2;
	end
end
yc = Yc(:, v)*Tc;

% matrices
E = .5*(yc'*Aplus*y + yc'*Aminus*y);
A = -.5/dCoef*(yc'*Aplus*y - yc'*Aminus*y);
if exist('b') == 1,
	B = yc'*b;
end
if exist('c') == 1,
	C = c*y;
end
