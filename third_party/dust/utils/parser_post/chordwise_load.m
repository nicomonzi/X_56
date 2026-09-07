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
% This matlab function parse the chordwise loads.
% The input is the same as the one required in the filed "basename" in the
% dust_post.in file. 
% The output is a data structure of the form:
%   chord_data_station.<field>.time     => simulation time
%   chord_data_station.<field>.value    => value matrix n_time x m_sec
%   chord_data_station.x_ref            => position of the center in chordwise direction
%                                          at first time step 
%   chord_data_station.z_ref            => position of the center in flapwise direction
%                                          at first time step 
% For further details regarding the <field> please see the input manual.

function chord_data = chordwise_load(file_out,n_station)

    field = {'Pres';
        'Cp';
        'dFx';
        'dFz';
        'dNx';
        'dNz';
        'dTx';
        'dTz';
        'x_cen';
        'z_cen'};

    % name file
    chord_data_station = struct('n_time',[],'n_chord',[],'x_ref',[],'z_ref',[],'chord_length',[], ...
        'spanwise_location',[], 'Pres',[],'Cp',[],'dFx',[], 'dFz',[],'dNx',[],'dNz',[], ...
        'dTx',[],'dTz',[],'x_cen',[],'z_cen',[]);
    
    for j = 1:n_station
        for i = 1:numel(field)
            file = [file_out '_' num2str(j) '_' field{i} '.dat'];
            chord_data_station = parse_file_chordwise(file, field{i}, chord_data_station);
        end
        chord_data(j) = chord_data_station; 
    end

end


function chord_data_station = parse_file_chordwise(file, field, chord_data_station)

    fid = fopen(file, 'rt');

    % initialization
    it = 0;
    skip_line = false;
    chord_data_station.n_time = 0;
    while skip_line || ~feof(fid)

        if ~skip_line
            line_new = fgets(fid);
        else
            skip_line = false;
        end

        it = it + 1;

        if it == 2
            line_2 = textscan(line_new, '# spanwise_location: %5.3f ; chord_length: %5.3f');
            chord_data_station.spanwise_location = line_2{1};
            chord_data_station.chord_length = line_2{2};
        end

        if it == 3
            line_3 = textscan(line_new,'# n_chord : %d ; n_time : %d. Next lines: x_chord , z_chord');
            chord_data_station.n_chord = line_3{1};
            chord_data_station.n_time = line_3{2};
        end

        if it == 4
            line_4 = textscan(line_new,repmat('%n',[1, chord_data_station.n_chord]));
            chord_data_station.x_ref = cell2mat(line_4(1:end));
        end

        if it == 5
            line_5 = textscan(line_new,repmat('%n',[1, chord_data_station.n_chord]));
            chord_data_station.z_ref = cell2mat(line_5(1:end));
        end

        if it > 6
            while skip_line || ~feof(fid)

                if ~skip_line
                    line = fgets(fid);
                else
                    skip_line = false;
                end

                line = textscan(line, repmat('%n',[1, chord_data_station.n_chord + 1]));
                chord_data_station.(field).time(it-6,1) = line{1};
                chord_data_station.(field).value(it-6,:) = cell2mat(line(2:end));
                it = it + 1;

            end
        elseif it > 6 + chord_data_station.n_time
            break
        end
    end
end





