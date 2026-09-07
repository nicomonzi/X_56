%./\\\\\\\\\\\...../\\\......./\\\..../\\\\\\\\\..../\\\\\\\\\\\\\.
%.\/\\\///////\\\..\/\\\......\/\\\../\\\///////\\\.\//////\\\////..
%..\/\\\.....\//\\\.\/\\\......\/\\\.\//\\\....\///.......\/\\\......
%...\/\\\......\/\\\.\/\\\......\/\\\..\////\\.............\/\\\......
%....\/\\\......\/\\\.\/\\\......\/\\\.....\///\\...........\/\\\......
%.....\/\\\......\/\\\.\/\\\......\/\\\.......\///\\\........\/\\\......
%......\/\\\....../\\\..\//\\\...../\\\../\\\....\//\\\.......\/\\\......
%.......\/\\\\\\\\\\\/....\///\\\\\\\\/..\///\\\\\\\\\/........\/\\\......
%........\///////////........\////////......\/////////..........\///.......
%%=========================================================================
% 
% Copyright (C) 2018 - 2022 Politecnico di Milano,
% with support from A^3 from Airbus
% and Davide Montagnani,
% Matteo Tugnoli,
% Federico Fonte
% 
% This file is part of DUST, an aerodynamic solver for complex
% configurations.
% 
% Permission is hereby granted, free of charge, to any person
% obtaining a copy of this software and associated documentation
% files (the "Software"), to deal in the Software without
% restriction, including without limitation the rights to use,
% copy, modify, merge, publish, distribute, sublicense, and / or sell
% copies of the Software, and to permit persons to whom the
% Software is furnished to do so, subject to the following
% conditions:
% 
% The above copyright notice and this permission notice shall be
% included in all copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
% EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
% OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
% NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
% HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
% WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
% FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
% OTHER DEALINGS IN THE SOFTWARE.
% 
% Author: Alessandro Cocco, Alberto Savino
%=========================================================================
%
% This function transfers the rotor sectional loads of DUST to the
% SU2_VARIABLE_LOAD Disk model.
% The data structure is: 
% data.file_dust; % root of filename of dust sectional load
% data.rho;       
% data.rpm; 
% data.n_blade;
% data.hub; % [0. 0. 0.] % hub center 
% data.axis; % [1.0 0. 0.] axis direction 
% data.file_su2; % filename of SU2 output .dat file

function dust2su2_actdisk(data)

    file_out = data.file_dust; 
    rho = data.rho;
    rpm = data.rpm; 
    n_blade = data.n_blade;
    hub =  data.hub;    % [0. 0. 0.] % hub center 
    axis = data.axis;   % [1.0 0. 0.] axis direction 
    file_su2  = data.file_su2; 


    short_data = true; 
    sec_data = sec_load(file_out, short_data);
    radius = sec_data.sec(end) + sec_data.l_sec(end)/2; 
    r = sec_data.sec/radius;  
    dFz = mean(sec_data.Fz.value);
    dFx = mean(sec_data.Fx.value);
    
    den  = (2*rho*(rpm/60)^2*(radius*2)^2); 
    dt_dr = n_blade.*(dFz*pi.*r)./den;
    dq_dr = -n_blade.*(dFx*pi^2.*r.^2)./den; 
    mat2write = [r',dt_dr',dq_dr', zeros(numel(r),1)];
    
    %% write su2.dat file
    fileID = fopen(file_su2,'wt');
    fprintf(fileID, '#\n#\n#\n#\n#\n#\n#\n#\n'); 
    fprintf(fileID, 'MARKER_ACTDISK= DISK DISK_OUT\n'); 
    fprintf(fileID, 'CENTER= %f %f %f\n', hub(1), hub(2), hub(3)); 
    fprintf(fileID, 'AXIS= %f %f %f\n', axis(1), axis(2), axis(3)); 
    fprintf(fileID, 'RADIUS= %f\n', radius); 
    fprintf(fileID, 'RPM= %f\n', rpm); 
    fprintf(fileID, 'NROW= %f\n', numel(r)); 
    fprintf(fileID, '# rs=r/R    dCT/drs     dCP/drs     dCR/drs\n'); 
    for i = 1:numel(r)
        fprintf(fileID, '%10.4f%10.4f%10.4f%10.4f\n', mat2write(i,:)); 
    end
    
    fclose(fileID);
    
    %% plot 
    subplot(1,2,1)
    plot(r, dt_dr)
    title('dT/dr')
    subplot(1,2,2)
    plot(r,dq_dr)
    title('dQ/dr')

end
