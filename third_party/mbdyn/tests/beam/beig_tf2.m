beig4_000;

% fmax = 400;
fmax = 3000;

b = zeros(127, 1);
b(12*4 + 9, 1) = -1;

c = zeros(1, 127);
c(1, 12*4 + 3) = 1;

vv = [100  101  106  107];
vvc = [87  88  93  94];

mbeig;

dw = 5.;
for i = 1:1000,
	w(i) = (i - 1)*dw;
	xx(:, i) = (j*w(i)*E - A)\B;
	yy(:, i) = C*xx(:, i);
end

% residualization...
if 0,
	[U, S, V] = svd(A\E);

	vd = [1 2];
	vs = [3 4];

	Sd = S(vd, vd);

	Ud = U(:, vd);
	Us = U(:, vs);
	Vd = V(:, vd);
	Vs = V(:, vs);

	Add = Ud'*Vd;
	Ads = Ud'*Vs;
	Asd = Us'*Vd;
	Ass = Us'*Vs;

	Bd = Ud'*(A\B);
	Bs = Us'*(A\B);

else
	[U, S, V] = svd(E);

	vd = [3 4];
	vs = [1 2];

	Sd = S(vd, vd);

	Ud = U(:, vd);
	Us = U(:, vs);
	Vd = V(:, vd);
	Vs = V(:, vs);

	Add = Ud'*A*Vd;
	Ads = Ud'*A*Vs;
	Asd = Us'*A*Vd;
	Ass = Us'*A*Vs;

	Bd = Ud'*B;
	Bs = Us'*B;
end

Cd = C*Vd;
Cs = C*Vs;

AA = Sd\(Add - Ads*(Ass\Asd));
BB = Sd\(Bd - Ads*(Ass\Bs));
CC = Cd - Cs*(Ass\Asd);
DD = -Cs*(Ass\Bs);

dw = 5.;
for i = 1:1000,
	w(i) = (i - 1)*dw;
	xi(:, i) = (j*w(i)*eye(size(AA)) - AA)\BB;
	yi(:, i) = CC*xi(:, i) + DD;
end

