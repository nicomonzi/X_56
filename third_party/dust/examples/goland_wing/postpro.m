% FLUTTER Goland wing
clear; close all; clc; 

% nc file name and path
file_nc = fullfile('structure','Output.nc');

% forces and moments in the first clamped node (wing root)
force  = ncread(file_nc, 'elem.joint.0.F');
moment =  ncread(file_nc, 'elem.joint.0.M');
time   =  ncread(file_nc, 'time');

% displacement and rotation of the last tip node
node_tip(:,:) =  ncread(file_nc, 'node.struct.8.X');
rot_tip  =  ncread(file_nc, 'node.struct.8.Phi');

% plot
figure()
plot(time, node_tip(3,:), 'LineWidth', 2)
xlabel('time [sec]')
ylabel('tip vertical displacement [m]')

