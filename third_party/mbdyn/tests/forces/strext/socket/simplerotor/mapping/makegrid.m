%
% MBDyn (C) is a multibody analysis code. 
% http://www.mbdyn.org
% 
% Copyright (C) 1996-2017
% 
% Pierangelo Masarati   <masarati@aero.polimi.it>
% Paolo Mantegazza      <mantegazza@aero.polimi.it>
% 
% Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
% via La Masa, 34 - 20156 Milano, Italy
% http://www.aero.polimi.it
% 
% Changing this copyright notice is forbidden.
% 
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation (version 2 of the License).
% 
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
% 
% Author: Giuseppe Quaranta <quaranta@aero.polimi.it>
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Creates a plane grid for the example simplerotor.
% The grid represents the mean surface of a straight rotor blade composed by
% symmetric airfoils, with constant chord and constant twist.
% The twist is applied around the 1/4-chord axis. 	 	
% 
% Input:
%	r0  	cutout radius. 		DEFAULT 0.25
%	rf  	external blade radius.	DEFAULT 5.0
%	chord   blade chord.		DEFAULT 0.5
%	sv	blade twist in deg.	DEFAULT 5.0
%	nc      number of elements in chord.	DEFAULT 10
%	ns	number of elements in span.	DEFAULT 20
%
% The grid data are saved in the file basicgrid.dat	
function makegrid(varargin)

r0 = 0.25;
rf = 5.0;
chord = 0.5;
sv = 5;
nc = 10;
ns = 20;

if nargin > 0 && ~isempty(varargin{1})
	r0 = varargin{1};
end
if nargin > 1 && ~isempty(varargin{2})
	rf = varargin{2};
end
if nargin > 2 && ~isempty(varargin{3})
	chord = varargin{3};
end
if nargin > 3 && ~isempty(varargin{4})
	sv = varargin{4};
end
if nargin > 4 && ~isempty(varargin{5})
	nc = varargin{5};
end
if nargin > 5 && ~isempty(varargin{6})
	ns = varargin{6};
end

sv = pi/180*sv;
t_nc = 2;
t_ns = 2;

dr = rf - r0;
npt = nc * ns;
ngt = (nc +1) * (ns + 1);

pos = zeros(3,ngt);
conn = zeros(npt,4);

switch t_nc
    case 1
        % suddivsione costante
        sudd_corda = 0:1/nc:1.;
    case 2
        % suddvisione sinusoidale solo Ba
        sudd_corda = 1 - cos(pi/2.*(0:1/nc:1));
    case 3
        % suddvisione sinusoidale ba e bu
        sudd_corda = 1/2*(1 - cos(pi.*(0:1/nc:1)));
    case 4 
        sudd_corda = 1/(10*nc) : (1-1/10)/(nc^2-1) : 1/nc;
        sudd_corda = [0.0, sudd_corda /sum(sudd_corda)];
        for i = 2:length(sudd_corda)
            sudd_corda(i) = sudd_corda(i) + sudd_corda(i-1);
        end
end

switch t_ns
    case 1
        % suddivsione costante
        sudd_span = 0:1/ns:1.;
    case 2
        % suddvisione sinusoidale solo tip
        sudd_span = sin(pi/2.*(0:1/ns:1));
    case 3
        % suddvisione sinusoidale root e tip
        sudd_span = 1/2*(1 - cos(pi.*(0:1/ns:1)));
        
    case 4 
        sudd_span = 1/ns : (1/10-1)/(ns^2-1) : 1/(10*ns);
        sudd_span = [0.0, sudd_span /sum(sudd_span)];
        for i = 2:length(sudd_span)
            sudd_span(i) = sudd_span(i) + sudd_span(i-1);
        end
end


for i = 1:ns + 1
	R = [ cos(sudd_span(i) * sv), 0, -sin(sudd_span(i) * sv) ;
       		 0                  , 1,               0         ; 
      	      sin(sudd_span(i) * sv), 0,  cos(sudd_span(i) * sv)];
    pm = [0; r0 + dr*sudd_span(i); 0]; 
    ba = pm + R*[-.25*chord; 0; 0];
    bu = pm + R*[.75*chord; 0; 0];
    dc = (bu - ba);
    for j = 1: nc + 1
        pos(:,(i-1)*(nc +1)+j) = ba + sudd_corda(j) * dc;
    end
end

n01 = 0;
np = 1;
for i=1:ns
    for j = 1:nc
        n01 = n01 + 1;
        n02 = n01 + 1;
        n03 = n02 + (nc + 1);
        n04 = n03 - 1;
        conn(np,:) = [n01, n02, n03, n04];
        np = np + 1;
    end
    n01 = n01 + 1;
end
 
disp('Printing grid file...');
fid = fopen('basicgrid.dat','w');
fprintf(fid, '# Program NUVOLA GRID\n');
fprintf(fid, '%4d\n',ngt);
for i=1:ngt
    fprintf(fid, '%16.10f%16.10f%16.10f\n', pos(1,i), pos(2,i), pos(3,i));
end
fprintf(fid, '%4d\n',npt);
for i=1:npt
    fprintf(fid, '%8d%8d%8d%8d\n', conn(i,1), conn(i,2), conn(i,3), conn(i,4));
end

fclose(fid);
disp('...completed');
