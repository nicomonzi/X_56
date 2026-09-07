% $Header: /var/cvs/mbdyn/mbdyn/mbdyn-1.0/tests/genel/hopf-bifurcation/hb.m,v 1.2 2012/04/19 09:28:17 masarati Exp $
%
% integrate the Hopf bifurcation problem

opt = odeset("RelTol", 1e-3, ...
 	"AbsTol", 1e-6, ...
 	"InitialStep", 1e-2, ...
 	"MaxStep", 1e-2);

data.lambda = .1;
data.alpha = -.01;
data.beta = 1.;
data.b = data.alpha + sqrt(-1)*data.beta;

data.r = sqrt(-data.lambda/data.alpha);
data.omega = data.beta*data.r^2;

% x0 = [4; 0];
x0 = [1; 0];

sol = ode45(@hbf, [0; 100], x0, opt, data);

