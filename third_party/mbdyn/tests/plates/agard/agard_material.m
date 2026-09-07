%proprieta' legno
%E_l = 3.151E9;
%E_t = 4.162E8;
%G = 4.392E8;
%nu_lt = 0.31;
syms E_l E_t G nu_lt positive
nu_tl = nu_lt / E_l * E_t;

%coefficienti per taglio trasverso e per rigidezza drilling
%as = 1.;
%at = 0.01;
syms as at positive

%rigidezza materiale stato piano sforzo
D = 1 /( 1 - nu_lt * nu_tl) * [ ...
	E_l E_l*nu_tl 0; ... 
	E_t * nu_lt E_t 0; ... 
	0 0 G * ( 1 - nu_lt * nu_tl) ...
];

%spessore legno
% t = 0.1;
syms t positive

%rigidezza da CLT
Ewoodblock = [ ...
	t*D zeros(3); ...
	zeros(3) 1/12*t^3*D ...
];

format short g

Ewoodblock

% Eshell = zeros(12, 12);

Incidence_Eshell = [1 2 4 5 7 8 10 11];

Perm_Ewoodblock = [1 3 3 2 4 6 6 5]; 

Shear_multiplier = diag([1 2 1 2 1 1 1 2 1 2 1 1]);

Eshell(Incidence_Eshell , Incidence_Eshell) = Ewoodblock(Perm_Ewoodblock, Perm_Ewoodblock);

Eshell(3, 3) = Eshell(2, 2) * as;
Eshell(6, 6) = Eshell(2, 2) * as;
Eshell(9, 9) = Eshell(8, 8) * at;
Eshell(12, 12) = Eshell(8, 8) * at;

Eshell = Eshell * Shear_multiplier
