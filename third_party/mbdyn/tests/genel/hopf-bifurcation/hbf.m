% $Header: /var/cvs/mbdyn/mbdyn/mbdyn-1.0/tests/genel/hopf-bifurcation/hbf.m,v 1.1 2009/07/06 20:44:42 masarati Exp $
%
% integrate the Hopf bifurcation problem

function ydot = hbf(tvar, yvar, varargin)

data = varargin{1};

z2 = yvar'*yvar;

ydot(1, 1) = yvar(1)*(data.lambda + data.alpha*z2) - yvar(2)*data.beta*z2;
ydot(2, 1) = yvar(1)*data.beta*z2 + yvar(2)*(data.lambda + data.alpha*z2);

endfunction
