clear all, close all, clc

%hinge moment check
%MBDyn nc
file_nc = fullfile('coupled','structure','Output','test.nc');
t_mb = ncread(file_nc, 'time');
hm_mb_r = ncread(file_nc, 'elem.joint.201.m');
hm_mb_l = ncread(file_nc, 'elem.joint.301.m');
% DUST postpro
file_du_r = fullfile('dust','Postpro','test_hm_Wing_Flap_r.dat');
file_du_r = dlmread(file_du_r, '' , 4, 0);
file_du_l = fullfile('dust','Postpro','test_hm_Wing_Flap_l.dat');
file_du_l = dlmread(file_du_l, '' , 4, 0);

figure()
set(gcf,'Position',[100 100 850 450])
subplot(1,2,1)
plot(t_mb, hm_mb_r(2,:), 'LineWidth', 2)
hold on, grid on
plot(file_du_r(:,1), file_du_r(:,6), 'LineWidth', 2)
title('right flap')
legend('nc file','DUST')
xlabel('time [sec]')
ylabel('h_m [Nm]')

subplot(1,2,2)
plot(t_mb, hm_mb_l(2,:), 'LineWidth', 2)
hold on, grid on
plot(file_du_l(:,1), file_du_l(:,6), 'LineWidth', 2)
title('left flap')
legend('nc file','DUST')
xlabel('time [sec]')
ylabel('h_m [Nm]')