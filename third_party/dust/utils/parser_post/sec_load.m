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
% Author: Alessandro Cocco 
%=========================================================================
%
% This matlab function parse the sectional loads. 
% The input is the same as the one required in the filed "basename" in the 
% dust_post.in file. 
% sec_data = sec_load(file_out, short_data)
% file_out: fullpath of the output.dat files without extension and attributes  
% If short_data is true,  then the parsed fields are Fx, Fy, Fz, Mo. 
% The output is a data structure of the form: 
%   sec_data.<field>.time     => simulation time 
%                   .value    => value matrix n_time x m_sec 
%                   .X        => position from global to local 
%                   .R        => rotation from global to local 
%   sec_data.sec              => position of the section from which the 
%                                <field> are extracted 
% For further details regarding the <field> please see the input manual. 

function sec_data = sec_load(file_out, short_data)
        
    if short_data
        field = {'Mo';
        'Fz';
        'Fy';
        'Fx';
        };
    else
        field = {'vel_outplane';
        'vel_2d';
        'up_x';
        'up_y';
        'up_z';
        'Mo';
        'Fz';
        'Fy';
        'Fx';
        'Cm';
        'Cl';
        'Cd';
        'alpha'};
    end
        
    

    % name file
    sec_data = struct('n_time', [], 'n_sec', [], 'sec', [], 'l_sec', [], 'chord', [], 'vel_outplane', [] ...
        , 'vel_outplane_isolated', [], 'vel_2d_isolated', [], 'vel_2d', [], 'up_x', [], 'up_y', [], 'up_z', [], ...
        'Mo',[],'Fx',[],'Fy',[],'Fz',[],'Cm',[],'Cl',[],'Cd',[],'alpha_isolated',[],'alpha',[]);

    for i = 1:numel(field)
        file = [file_out '_' field{i} '.dat'];
        sec_data = parse_file_sectional(file, field{i}, sec_data);
    end

end


function sec_data = parse_file_sectional(file, field, sec_data)

    fid = fopen(file, 'rt');
    
    % initialization
    it = 0;
    skip_line = false;
    sec_data.n_time = 0; 
    long_field = {'Mo';
                'Fy';
                'Fz';
                'Fx'};
    while skip_line || ~feof(fid)
    
        if ~skip_line
            line_new = fgets(fid);
        else
            skip_line = false;
        end
    
        it = it + 1;
    
        if it == 2
            line_2 = textscan(line_new, '# n_sec : %d ; n_time : %d. Next lines: y_cen , y_span, chord');
            sec_data.n_sec = line_2{1};
            sec_data.n_time = line_2{2};
        end
    
        if it == 3
            line_3 = textscan(line_new,[repmat('%n',[1, sec_data.n_sec])]);
            sec_data.sec = cell2mat(line_3(1:end));
        end
    
        if it == 4
            line_4 = textscan(line_new,[repmat('%n',[1, sec_data.n_sec])]);
            sec_data.l_sec = cell2mat(line_4(1:end));
        end

        if it == 5
            line_5 = textscan(line_new,[repmat('%n',[1, sec_data.n_sec])]);
            sec_data.chord = cell2mat(line_5(1:end));
        end

    
        if it > 6
            while skip_line || ~feof(fid)
    
                if ~skip_line
                    line = fgets(fid);
                else
                    skip_line = false;
                end
    
                if any(strcmpi(field,long_field))
                    line = textscan(line, repmat('%n',[1, 1 + sec_data.n_sec + 9 + 3]));
                    sec_data.(field).time(it-6,1) = line{1};
                    sec_data.(field).value(it-6,:) = cell2mat(line(2:end-12));
                    sec_data.(field).R(:,:,it-6) = reshape(cell2mat(line(end-11:end-3)),3,3);
                    sec_data.(field).X(:,:,it-6) = line(end-2:end);
                    it = it + 1;
                else
                    line = textscan(line, repmat('%n',[1, sec_data.n_sec + 1]));
                    sec_data.(field).time(it-6,1) = line{1};
                    sec_data.(field).value(it-6,:) = cell2mat(line(2:end));
                    it = it + 1;
                end
    
            end
        elseif it > 6 + sec_data.n_time
            break
        end
    end
end





