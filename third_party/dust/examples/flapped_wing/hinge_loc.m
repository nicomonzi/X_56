clear all, close all, clc

sweep = 30;
twist = 0;
dihed = 10;
n1 = [0.5 2 0]';
n2 = [0.5 4 0]';

R =@(x,y,z) [cosd(z) -sind(z) 0; sind(z) cosd(z) 0; 0 0 1]*[cosd(y) 0 sind(y); 0 1 0;-sind(y) 0 cosd(y)]*[1 0 0; 0 cosd(x) -sind(x); 0 sind(x) cosd(x)];
rot = R(dihed,twist,-sweep);

Node1 = rot*n1
Node2 = rot*n2