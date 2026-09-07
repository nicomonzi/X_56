clear; close all; clc; 
rr = dlmread(fullfile('Output','test_mesh_points_0000.dat'));
ee = dlmread(fullfile('Output','test_mesh_elems_0000.dat'));

%plot3(rr(:,1),rr(:,2),rr(:,3), '-*')
patch('Faces',ee,'Vertices',rr,'FaceColor',[1 1 1])
axis equal