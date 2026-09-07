% $Header: /var/cvs/mbdyn/mbdyn/mbdyn-1.0/tests/beam/DETC2007-35511/DETC2007_35511.m,v 1.1 2009/08/21 12:04:29 masarati Exp $
%
% Copyright 2009 Pierangelo Masarati <masarati@aero.polimi.it>
%
% Space manipulator based on DETC2007-35511
% Proceedings of IDETC/CIE 2007
% ASME 2007 International Design Engineering Technical Conferences &
% Computers and Information in Engineering Conference
% September 4-7, 2007, Las Vegas, USA
% INVESTIGATION OF BOUNDARY CONDITIONS FOR FLEXIBLE MULTIBODY
% SPACECRAFT DYNAMICS (DRAFT VERSION)
% John R. MacLean
% An Huynh
% Leslie J. Quiocho
%
% POSSIBLE INCONSISTENCIES
%
% 1) XI[4] and OMEGA[4] in Table 2 do not conform to
%	OMEGA = (MU/L)^2*sqrt(EI/RHO)
%
% 2) torque levels from Figure 7 (+/- ~44 Nm) looks overkill; divided by 9.81
%    yields comparable results

RHO = 2.82;		% kg/m
EI = 5.974^2*RHO;	% Nm^2
L = 3.372;		% m
XI = 5.18;		% [adim]
M_BEAM = RHO*L;		% kg
M_TIP = XI*M_BEAM;	% kg

m_1 = M_BEAM + M_TIP;
X_CM = (M_BEAM*L/2 + M_TIP*L)/m_1;
c_1_cross = [0., 0., 0.;
	     0., 0., -X_CM;
	     0., X_CM, 0.];
J_1 = M_BEAM*L^2/12 + M_BEAM*(X_CM - L/2)^2 + M_TIP*(L - X_CM)^2;

M_rr_1 = [m_1, m_1*X_CM;
	  m_1*X_CM, J_1];

% NOTE: the paper reports
% OMEGA = 5.28        20.9        46.7        83.4
% MU =    3.171011466 6.298221427 9.434866134      15.71404495
% there is clearly a mismatch, as MU=15.71404495 is not related to OMEGA=83.4
% provisionally fixing...
OMEGA = [5.28; 20.9; 46.7; 83.4; 129.73766239713339];
MU = [3.171011466; 6.298221427; 9.434866134; 12.59905529960816; 15.71404495];

SIGMA = EI/RHO*(MU/L).^4;

A0 = sqrt(3/(L^3*(1 + 3*XI)));
A = 1./sqrt(L*(1 + XI + 2*XI*MU.*(XI*MU + coth(MU))));

Int_X = M_BEAM*(A./MU).*((cosh(MU) - 1)./sinh(MU) - (cos(MU) - 1)./sin(MU));
Int_rX = (A./sinh(MU)).*(L^2*cosh(MU)./MU - L^2*sinh(MU)./(MU.^2)) ...
	+ (A./sin(MU)).*(-L^2*cos(MU)./MU + L^2*sin(MU)./(MU.^2));

M_re = [Int_X'; Int_rX'];

M_ee = -(M_BEAM/L)*4*XI*A*A';
M_ee = M_ee - diag(diag(M_ee));
M_ee = M_ee + diag((M_BEAM/L)*((1 - 3*XI + 2*XI*MU.*(XI*MU + coth(MU)))./(1 + XI + 2*XI*MU.*(XI*MU + coth(MU)))));

K_ee = EI*diag((MU/L).^4);

x = [0:.01:1]*L;

X = (A*ones(size(x))).*(sinh((MU/L)*x)./sinh(MU*ones(size(x))) + sin((MU/L)*x)./sin(MU*ones(size(x))));
ALPHA = ((A.*MU/L)*ones(size(x))).*(cosh((MU/L)*x)./sinh(MU*ones(size(x))) + cos((MU/L)*x)./sin(MU*ones(size(x))));

X_0 = X(:, 1);
X_L = X(:, end);
ALPHA_0 = ALPHA(:, 1);
ALPHA_L = ALPHA(:, end);


