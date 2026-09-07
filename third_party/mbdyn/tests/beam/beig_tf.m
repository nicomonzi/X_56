beig4_000;

% fmax = 400;
fmax = 3000;

b = zeros(127, 1);
b(12*4 + 9, 1) = 1;

c = zeros(1, 127);
c(1, 12*4 + 3) = 1;

mbeig;

dw = 5.;
for i = 1:1000,
	w(i) = (i - 1)*dw;
	xx(:, i) = (j*w(i)*E - A)\B;
	yy(:, i) = C*xx(:, i);
end

% residualization...
[U, S, V] = svd(E);
Sd = S(1:2, 1:2);
Ud = U(:, 1:2);
Us = U(:, 3:4);
Vd = V(:, 1:2);
Vs = V(:, 3:4);

Add = Ud'*A*Vd;
Ads = Ud'*A*Vs;
Asd = Us'*A*Vd;
Ass = Us'*A*Vs;

Bd = Ud'*B;
Bs = Us'*B;

Cd = C*Vd;
Cs = C*Vs;

AA = Sd\(Add - Ads*(Ass\Asd));
BB = Sd\(Bd - Ads*(Ass\Bs));
CC = Cd - Cs*(Ass\Asd);
DD = -Cs*(Ass\Bs);

