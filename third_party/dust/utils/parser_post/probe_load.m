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
% This matlab function parse the probe fields.
% The input is the same as the one required in the filed "basename" in the
% dust_post.in file.
% The field can be velocity, pressure, cp
% The output is a data structure of the form:
%   probe_data.<field>.time     => simulation time
%                     .value    => value matrix field_dim x n_time x
%                     n_probe
%   probe_data.n_probe          => number of probe points 
%   probe_data.point            => coordinates matrix 3 x n_probe

function probe_data = probe_load(file)

    fid = fopen(file, 'rt');
    
    % initialization
    it = 0;
    it_mat = 0; 
    var = 0;
    skip_line = false;
    velocity = false;
    pressure = false;
    cp = false;
    probe_data.n_time = 0;
    while skip_line || ~feof(fid)
    
        if ~skip_line
            line_new = fgets(fid);
        else
            skip_line = false;
        end
    
        it = it + 1;
            
        if it == 1
            line_1 = textscan(line_new,' # N. of point probes:           %d');
            probe_data.n_probe = line_1{1};
        end

        if it == 2
            line_2 = textscan(line_new,[repmat('%n',[1, probe_data.n_probe])]);
            probe_data.point(1,:) = cell2mat(line_2);
        end
    
        if it == 3
            line_3 = textscan(line_new,[repmat('%n',[1, probe_data.n_probe])]);
            probe_data.point(2,:) = cell2mat(line_3);
        end
    
        if it == 4
            line_4 = textscan(line_new,[repmat('%n',[1, probe_data.n_probe])]);
            probe_data.point(3,:) = cell2mat(line_4);
        end
    
        if it == 5
            line_5 = textscan(line_new, '# n_time: %d');
            probe_data.n_time = line_5{1};
        end
    
        if it == 6
            line_6 = textscan(line_new,'#    t     %d (  ux  uy  uz  p  cp )');
            probe_data.n_probe = line_6{1};
            if contains(line_new,'ux')
                var = var + 3;
                velocity = true;
                probe_data.velocity.time = zeros(probe_data.n_time,1); 
                probe_data.velocity.value = zeros(3,probe_data.n_time, probe_data.n_probe); 
            end
            if contains(line_new,'p')
                var = var + 1;
                pressure = true;
                probe_data.pressure.time = zeros(probe_data.n_time,1); 
                probe_data.pressure.value = zeros(1,probe_data.n_time, probe_data.n_probe); 
            end
            if contains(line_new,'cp')
                var = var + 1;
                cp = true;       
                probe_data.cp.time = zeros(probe_data.n_time, 1);
                probe_data.cp.value = zeros(1,probe_data.n_time, probe_data.n_probe);
            end

        end
    
    
        if it > 5
            while skip_line || ~feof(fid)
    
                if ~skip_line
                    line = fgets(fid);
                else
                    skip_line = false;
                end

                it_mat = it_mat + 1; 
                
                line = textscan(line, repmat('%n',[1, 1 + probe_data.n_probe*var]));
                matrix = reshape(cell2mat(line(2:end)), [var, probe_data.n_probe]);
                if var == 1
                    if pressure
                    probe_data.pressure.time(it_mat,1) =  line{1};
                        probe_data.pressure.value(1,it_mat,:) = ...
                            reshape(matrix(1,:),[1,1,probe_data.n_probe]);
    
                    elseif cp
                    probe_data.cp.time(it_mat,1) =  line{1};
                    probe_data.cp.value(1,it_mat,:) = ...
                            reshape(matrix(1,:),[1,1,probe_data.n_probe]);
    
                    end
                end
                if var == 2
                    probe_data.pressure.time(it_mat,1) =  line{1};
                        probe_data.pressure.value(1,it_mat,:) = ...
                            reshape(matrix(1,:),[1,1,probe_data.n_probe]);
                    probe_data.cp.time(it_mat,1) =  line{1};
                    probe_data.cp.value(1,it_mat,:) = ...
                            reshape(matrix(2,:),[1,1,probe_data.n_probe]);
                end
                if var == 3 % only velocity
                    if velocity
                    probe_data.velocity.time(it_mat,1) =  line{1};
                        probe_data.velocity.value(1:3,it_mat,:) = ...
                            reshape(matrix(1:3,:),[3,1,probe_data.n_probe]);
                    end
                elseif var == 4
                    
                    if velocity
                        probe_data.velocity.time(it_mat,1) =  line{1};
                        probe_data.velocity.value(1:3,it_mat,:) = ...
                            reshape(matrix(1:3,:),[3,1,probe_data.n_probe]);
                    end
                    if pressure
                        probe_data.pressure.time(it_mat,1) =  line{1};
                        probe_data.pressure.value(1,it_mat,:) = ...
                            reshape(matrix(4,:),[1,1,probe_data.n_probe]);
                    end
                    if cp
                        probe_data.cp.time(it_mat,1) =  line{1};
                        probe_data.cp.value(1,it_mat,:) = ...
                            reshape(matrix(4,:),[1,1,probe_data.n_probe]);
                    end
                else % all output
                    probe_data.velocity.time(it_mat,1) =  line{1};
                        probe_data.velocity.value(1:3,it_mat,:) = ...
                            reshape(matrix(1:3,:),[3,1,probe_data.n_probe]);
                    probe_data.pressure.time(it_mat,1) =  line{1};
                        probe_data.pressure.value(1,it_mat,:) = ...
                            reshape(matrix(4,:),[1,1,probe_data.n_probe]);
                    probe_data.cp.time(it_mat,1) =  line{1};
                    probe_data.cp.value(1,it_mat,:) = ...
                            reshape(matrix(5,:),[1,1,probe_data.n_probe]);    
    
                end
    

            end
        elseif it > 5 + probe_data.n_time
            break
        end
    end
end





