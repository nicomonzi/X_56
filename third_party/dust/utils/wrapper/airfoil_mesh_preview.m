%./\\\\\\\\\\\...../\\\......./\\\..../\\\\\\\\\..../\\\\\\\\\\\\\.
%.\/\\\///////\\\..\/\\\......\/\\\../\\\///////\\\.\//////\\\////..
%..\/\\\.....\//\\\.\/\\\......\/\\\.\//\\\....\///.......\/\\\......
%...\/\\\......\/\\\.\/\\\......\/\\\..\////\\.............\/\\\......
%....\/\\\......\/\\\.\/\\\......\/\\\.....\///\\...........\/\\\......
%.....\/\\\......\/\\\.\/\\\......\/\\\.......\///\\\........\/\\\......
%......\/\\\....../\\\..\//\\\...../\\\../\\\....\//\\\.......\/\\\......
%.......\/\\\\\\\\\\\/....\///\\\\\\\\/..\///\\\\\\\\\/........\/\\\......
%........\///////////........\////////......\/////////..........\///.......
%=========================================================================
%
% Copyright (C) 2018-2023 Politecnico di Milano,
%                           with support from A^3 from Airbus
%                    and  Davide   Montagnani,
%                         Matteo   Tugnoli,
%                         Federico Fonte
%
% This file is part of DUST, an aerodynamic solver for complex
% configurations.
%
% Permission is hereby granted, free of charge, to any person
% obtaining a copy of this software and associated documentation
% files (the "Software"), to deal in the Software without
% restriction, including without limitation the rights to use,
% copy, modify, merge, publish, distribute, sublicense, and/or sell
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
%
%=========================================================================

% This script serves as preview of the airfoil mesh that will be created
% in dust_pre after having loaded a .dat profile in Selig format

clear; close all; clc;

%% Input

rr = readmatrix('naca4412.dat'); % file containing the profile coordinate

n_elem = 200; % number of elements to discretize

% discretization type (see DUST manual)
%mesh_type = 'uniform';
%mesh_type = 'cosineLE';
%mesh_type = 'cosineTE';
%mesh_type = 'cosine';
%-> geometric series discretization possibilities
%mesh_type = 'geoseriesLE';
%mesh_type = 'geoseriesTE';
r = 1/5;
mesh_type = 'geoseries';
r_le = 1/8;
r_te = 1/10;
%mesh_type = 'geoseries';
xh = 0.5; % hinge chordwise position (adimensional)
theta = deg2rad(-0); % hinge deflection (negative downwards)
offset = 0.03; % size of the blending region
r_le_aft = 1/10; % growth ratio at leading edge of the fixed part
r_te_aft = 1/8; % growth ratio at trailing edge of the fixed part
r_le_fore = 1/5; % growth ratio at leading edge of the hinged part
r_te_fore = 1/8; % growth ratio at trailing edge of the hinged part

% just for verication
airfoil_str = 'NACA4412';
M = str2double(airfoil_str(5));
P = str2double(airfoil_str(6));
SS = str2double(airfoil_str(7:8));
[x_anal, y_anal] = naca4digit(M, P, SS, 1, n_elem);

%% Interpolate dat file using a combination of Hermite Spline and a circle sector
%% some checks
% find minx and maxx
[rr_xmin, id_min] = min(rr(:, 1));
[rr_xmax, id_max] = max(rr(:, 1));
% check which is the lower side
if rr(2, 2) > 0
    rr_low2interp = rr(id_max:id_min, :);
    rr_up2interp = rr(id_min:end, :);
else
    rr_up2interp = rr(id_max:id_min, :);
    rr_low2interp = rr(id_min:end, :);
end

% split

% combine vectors
rr_comb = [rr_up2interp; rr_low2interp];
[rr_sort(:, 1), idsort] = sort(rr_comb(:, 1));
rr_sort(:, 2) = rr_comb(idsort, 2);
% remouve double elemnts
rr_unique = unique(rr_sort, 'rows');

yp_interp = interp1(rr_up2interp(:, 1), rr_up2interp(:, 2), rr_unique(:, 1));
yl_interp = interp1(rr_low2interp(:, 1), rr_low2interp(:, 2), rr_unique(:, 1));
rr_low = [rr_unique(:, 1), yp_interp];

rr_up = [rr_unique(:, 1), yl_interp];
rr_up = unique(rr_up, "rows");
rr_low = unique(rr_low, "rows");

% build middle line
% id_up = floor(size(rr,1)/2);
% rr_up = rr(1:id_up + 1,:); % assuming Seilig format
% rr_low = rr(id_up + 1:end,:);
rr_mean(:, 1) = rr_low(:, 1);
rr_mean(:, 2) = (rr_up(:, 2) + rr_low(:, 2)) / 2;
% get read of weird behaviuor
xmean = linspace(0, 1, 10);
ymean = interp1(rr_mean(:, 1), rr_mean(:, 2), xmean);
figure
plot(rr_mean(:, 1), rr_mean(:, 2))
hold on
plot(rr_up(:, 1), rr_up(:, 2))
plot(rr_low(:, 1), rr_low(:, 2))
plot(xmean, ymean)

rr_mean = [xmean' ymean'];

%%
% get m and p
[m, idx_m] = max(rr_mean(:, 2)); % y_mean line max thickness
p = rr_mean(idx_m, 1); % x_mean line max thickness
% the mean line is given by the expression
% y_m = m/p^2(2px - x^2)
% The expression is related NACA profile, but should be valid for all profile
% close to the leading edge.

% Get analytical derivative of the mean line in x = 0
m_line = 2 * m / p; % need to get the center of the circle -> linearized approach
tang_mean = (rr_mean(end - 1, 2) - rr_mean(end, 2)) / ...
    (rr_mean(end - 1, 1) - rr_mean(end, 1));

if abs(tang_mean) < 0.02
    tang_mean = 0;
end

xq_mean = linspace(0, 1, 20000);
% start iterative loop on m_line
maxiter = 10;

for i = 1:maxiter
    yq_mean = hermite_spline(rr_mean(:, 1), rr_mean(:, 2), xq_mean, m_line, tang_mean);
    [m, idx_m] = max(yq_mean); % y_mean line max thickness
    p = xq_mean(idx_m); % x_mean line max thickness
    m_line_new = 2 * m / p; % need to get the center of the circle -> linearized approach

    if abs(m_line_new - m_line) < 1e-6
        break
    end

    m_line = m_line_new;
end

% simple condition on m
if m_line < 0.03
    m_line = 0; % symmetric profile
end

% calculate thickness in a simple manner
t = max(rr_up(:, 2)) + min(rr_low(:, 1));
%t = 0.12;
% start iterative loop on coef in order to have the last point belonging to
% the origial curve
% radius of the circle at the leading edge
coef = 1.10; % nacaprofile
circ.radius = coef * t ^ 2;
maxiter = 100;
xup_int = linspace(0., 1, 20000);
yup_int = interp1(rr_up(:, 1), rr_up(:, 2), xup_int);

for i = 1:maxiter

    % get center of the circle with the following conditions:
    % 1. passing through [0 0]
    % 2. center lying on the linearized middle line
    % 3. radius given by the NACA formula
    circ.x_c = sqrt(circ.radius ^ 2 / (m_line ^ 2 + 1));
    circ.y_c = sqrt(circ.radius ^ 2 - circ.x_c ^ 2);

    % Define the region that is approximated by a circle arc
    % - get point at +/- 37.5 deg (angle that works the best)
    alpha = atan(circ.y_c / circ.x_c); % angle of m
    beta_up = (180 - 37.5) * pi / 180;
    beta_low = (180 + 37.5) * pi / 180;
    circ.x_plus = circ.radius * cos(beta_up + alpha) + circ.x_c;
    circ.y_plus = circ.radius * sin(beta_up + alpha) + circ.y_c;
    circ.x_minus = circ.radius * cos(beta_low + alpha) + circ.x_c;
    circ.y_minus = circ.radius * sin(beta_low + alpha) + circ.y_c;
    % check intersection
    % loop x_up profile
    for j = 1:numel(xup_int) - 1

        if circ.x_plus >= xup_int(j) && circ.x_plus <= xup_int(j + 1)
            ind_x = j;
        end

    end

    if circ.y_plus >= yup_int(ind_x) && circ.y_plus <= yup_int(ind_x + 1)
        break
    else
        radius_old = circ.radius;
        circ.radius = sqrt((xup_int(ind_x) - circ.x_c) ^ 2 + ...
            (yup_int(ind_x) - circ.y_c) ^ 2);

        if radius_old > circ.radius % naca profile
            circ.radius = radius_old;
            break
        end

    end

    disp(i)
end

% get the derivative for tangent at leading edge (to use in the
% splining process
circ.tan_plus = tan(atan(-cot(beta_up)) + alpha);
circ.tan_minus = tan(atan(-cot(beta_low)) + alpha);
circ.tan_le_up = tan(pi / 2 + alpha);
circ.tan_le_low = tan(-pi / 2 + alpha);

% leading edge circle sector
theta_start = beta_up + alpha;
theta_end = beta_low + alpha;
circ.theta = linspace(theta_start, theta_end, 300000);
circ.x = circ.radius .* cos(circ.theta) + circ.x_c;
circ.y = circ.radius .* sin(circ.theta) + circ.y_c;
plot(circ.x, circ.y)
% tangent at trailing edge: backward difference
tang_te_up = (rr_up(end - 1, 2) - rr_up(end, 2)) / ...
    (rr_up(end - 1, 1) - rr_up(end, 1));
tang_te_low = (rr_low(end - 1, 2) - rr_low(end, 2)) / ...
    (rr_low(end - 1, 1) - rr_low(end, 1));

% spline interpolation upper part
for i = 1:size(rr_up, 1)

    if rr_up(i, 1) > circ.x_plus
        idx = i;
        break
    end

end

point_up = [0 0;
            rr_up(idx:end, :)]; % point from ,dat file

% replace first point with circle point
point_up(1, :) = [circ.x_plus, circ.y_plus];

y_fake_pt_le = circ.tan_plus * (circ.x_plus * 0.9 - circ.x_plus) + circ.y_plus;
y_fake_pt_te = tang_te_up * (point_up(end, 1) * 1.1 - point_up(end, 1)) + point_up(end, 2);

point_up_updated = [circ.x_plus(1, 1) * 0.9 y_fake_pt_le;
                    point_up;
                    point_up(end, 1) * 1.1 y_fake_pt_te];

num_points = 1000; % Density of blue chain points between two red points
chain_points_up = catmull_rom_chain(point_up_updated, num_points);

%xq_up = linspace(circ.x_plus, 1, 5000); % from circle_up to TE
% yq_up = hermite_spline(point_up(:,1), point_up(:,2), xq_up, circ.tan_plus, tang_te_up);

% spline interpolation lower part
for i = 1:size(rr_low, 1)

    if rr_low(i, 1) > circ.x_plus
        idx = i;
        break
    end

end

point_low = [0 0;
             rr_low(idx:end, :)]; % point from ,dat file
point_low(1, :) = [circ.x_minus, circ.y_minus];

y_fake_pt_le = circ.tan_minus * (circ.x_minus * 0.9 - circ.x_minus) + circ.y_minus;
y_fake_pt_te = tang_te_low * (point_low(end, 1) * 1.1 - point_low(end, 1)) + point_low(end, 2);

point_low_updated = [circ.x_minus(1, 1) * 0.9 y_fake_pt_le;
                     point_low;
                     point_low(end, 1) * 1.1 y_fake_pt_te];

chain_points_low = catmull_rom_chain(point_low_updated, num_points);

%xq_low = linspace(circ.x_minus, 1, 5000);
%yq_low = hermite_spline(point_low(:,1), point_low(:,2), xq_low, circ.tan_minus, tang_te_low);

%% reorganize vectors to prepare mesh subdivision
[min_x, idx_circ_min] = min(circ.x); % get lowest x (in general is <= 0)

circ.x_up = circ.x(1:idx_circ_min);
circ.y_up = circ.y(1:idx_circ_min);
circ.x_low = circ.x(idx_circ_min:end);
circ.y_low = circ.y(idx_circ_min:end);
circ.x_up = fliplr(circ.x_up);
circ.y_up = fliplr(circ.y_up);

% theta_start = beta_up + alpha;
% theta_end =  pi - alpha;
% circ.theta = linspace(theta_start, theta_end, 150000);
% circ.x_up = circ.radius.*cos(circ.theta) + circ.x_c;
% circ.y_up = circ.radius.*sin(circ.theta) + circ.y_c;
% circ.x_up = fliplr(circ.x_up);
% circ.y_up = fliplr(circ.y_up);
%
% theta_end = beta_low + alpha;
% theta_start =  pi - alpha;
% circ.theta = linspace(theta_start, theta_end, 150000);
% circ.x_low = circ.radius.*cos(circ.theta) + circ.x_c;
% circ.y_low = circ.radius.*sin(circ.theta) + circ.y_c;

% upper
rr2mesh_up = [circ.x_up chain_points_up(2:end, 1)';
              circ.y_up chain_points_up(2:end, 2)'];
% get length and alpha
alpha_up = zeros(1, size(rr2mesh_up, 2));
dl_up = zeros(1, size(rr2mesh_up, 2));
alpha_up(1) = atan(circ.tan_le_up);
dl_up(1) = 0;

for i = 2:size(rr2mesh_up, 2)
    alpha_up(i) = atan((rr2mesh_up(2, i) - rr2mesh_up(2, i - 1)) / ...
        (rr2mesh_up(1, i) - rr2mesh_up(1, i - 1)));
    dl_up(i) = dl_up(i - 1) + sqrt((rr2mesh_up(1, i) - rr2mesh_up(1, i - 1)) ^ 2 + ...
        (rr2mesh_up(2, i) - rr2mesh_up(2, i - 1)) ^ 2);
end

% lower part
rr2mesh_low = [circ.x_low chain_points_low(2:end, 1)';
               circ.y_low chain_points_low(2:end, 2)'];

% get length and alpha
alpha_low = zeros(1, size(rr2mesh_low, 2));
dl_low = zeros(1, size(rr2mesh_low, 2));
alpha_low(1) = atan(circ.tan_le_low);
dl_low(1) = 0;

for i = 2:size(rr2mesh_low, 2)
    alpha_low(i) = atan((rr2mesh_low(2, i) - rr2mesh_low(2, i - 1)) / ...
        (rr2mesh_low(1, i) - rr2mesh_low(1, i - 1)));
    dl_low(i) = dl_low(i - 1) + sqrt((rr2mesh_low(1, i) - rr2mesh_low(1, i - 1)) ^ 2 + ...
        (rr2mesh_low(2, i) - rr2mesh_low(2, i - 1)) ^ 2);
end

% mean line
rr2mesh_mean = [xq_mean; yq_mean];
alpha_mean = zeros(1, size(rr2mesh_mean, 2));
dl_mean = zeros(1, size(rr2mesh_mean, 2));
alpha_mean(1) = m_line;
dl_mean(1) = 0;

for i = 2:size(rr2mesh_mean, 2)
    alpha_mean(i) = atan((rr2mesh_mean(2, i) - rr2mesh_mean(2, i - 1)) / ...
        (rr2mesh_mean(1, i) - rr2mesh_mean(1, i - 1)));
    dl_mean(i) = dl_mean(i - 1) + sqrt((rr2mesh_mean(1, i) - rr2mesh_mean(1, i - 1)) ^ 2 + ...
        (rr2mesh_mean(2, i) - rr2mesh_mean(2, i - 1)) ^ 2);
end

%% Mesh possibilites:
% uniform
% cosineLE
% cosineTE
% cosine
% geometric series
%   - LE
%   - TE
%   - both
%   - hinge

switch mesh_type
    case 'uniform'
        dcsi = linspace(0., 1, n_elem +1);
    case 'cosineLE'
        dcsi = linspace(0., 1, n_elem +1);
        dcsi = 1.0 - cos(0.5 * pi .* dcsi);
    case 'cosineTE'
        dcsi = linspace(0., 1, n_elem +1);
        dcsi = sin(0.5 * pi .* dcsi);
    case 'cosine'
        dcsi = linspace(0., 1, n_elem +1);
        dcsi = 1/2.0 .* (1.0 - cos(pi .* dcsi));
    case 'geoseriesLE'
        [~, dcsi] = geoseries(0., 1, n_elem, r);
    case 'geoseriesTE'
        [dcsi, ~] = geoseries(0., 1, n_elem, r);
    case 'geoseries'
        dcsi = geoseries_both(0., 1, n_elem, r_le, r_te);
    case 'geoseriesHI'
        dcsi = geoseries_hinge(0., 1, n_elem, xh, ...
            r_le_aft, r_te_aft, r_le_fore, r_te_fore);
end

dcsi_up = rescale_csi(dcsi, 0., dl_up(end));
x_up_csi = zeros(1, numel(dcsi_up));
y_up_csi = zeros(1, numel(dcsi_up));
x_up_csi(1) = rr2mesh_up(1, 1);
y_up_csi(1) = rr2mesh_up(2, 1);

for i = 1:numel(dcsi_up)

    for j = 1:numel(dl_up) - 1
        % two consecutive points x for which the point xq belongs
        if (dcsi_up(i) >= dl_up(j)) && (dcsi_up(i) <= dl_up(j + 1))
            x_up_csi(i) = rr2mesh_up(1, j) + cos(alpha_up(j)) * (dcsi_up(i) - dl_up(j));
            y_up_csi(i) = rr2mesh_up(2, j) + sin(alpha_up(j)) * (dcsi_up(i) - dl_up(j));
        end

    end

end

dcsi_low = rescale_csi(dcsi, 0., dl_low(end));
x_low_csi = zeros(1, numel(dcsi_low));
y_low_csi = zeros(1, numel(dcsi_low));
x_low_csi(1) = rr2mesh_low(1, 1);
y_low_csi(1) = rr2mesh_low(2, 1);

for i = 1:numel(dcsi_low)

    for j = 1:numel(dl_low) - 1
        % two consecutive points x for which the point xq belongs
        if (dcsi_low(i) >= dl_low(j)) && (dcsi_low(i) <= dl_low(j + 1))
            x_low_csi(i) = rr2mesh_low(1, j) + cos(alpha_low(j)) * (dcsi_low(i) - dl_low(j));
            y_low_csi(i) = rr2mesh_low(2, j) + sin(alpha_low(j)) * (dcsi_low(i) - dl_low(j));
        end

    end

end

% correct ends

dcsi_mean = rescale_csi(dcsi, 0., dl_mean(end));
x_mean_csi = zeros(1, numel(dcsi_mean));
y_mean_csi = zeros(1, numel(dcsi_mean));

for i = 1:numel(dcsi_mean)

    for j = 1:numel(dl_mean) - 1
        % two consecutive points x for which the point xq belongs
        if (dcsi_mean(i) >= dl_mean(j)) && (dcsi_mean(i) <= dl_mean(j + 1))
            x_mean_csi(i) = rr2mesh_mean(1, j) + cos(alpha_mean(j)) * (dl_mean(j + 1) - dl_mean(j));
            y_mean_csi(i) = rr2mesh_mean(2, j) + sin(alpha_mean(j)) * (dl_mean(j + 1) - dl_mean(j));
        end

    end

end

% linear interpolation to get y coordinates
% y_mesh_up = interp1(rr2mesh_up(1,:), rr2mesh_up(2,:), dcsi,'linear','extrap');
% y_mesh_low = interp1(rr2mesh_low(1,:), rr2mesh_low(2,:), dcsi,'linear','extrap');

% mean line
% y_mesh_mean = (y_mesh_low + y_mesh_up)/2;
rr_mean = [x_mean_csi; y_mean_csi];

rr_final = [fliplr(x_low_csi) x_up_csi(2:end);
            fliplr(y_low_csi) y_up_csi(2:end)];

% hinge deflection
rr_deflected = hinge_deflection(rr_final, theta, xh, offset);

%% plotting results
figure
plot(rr_final(1, :), rr_final(2, :), 'r-+', 'LineWidth', 2)
hold on
plot(rr(:, 1), rr(:, 2))
plot(rr_mean(1, :), rr_mean(2, :), 'r-+', 'LineWidth', 2)
plot(rr_deflected(1, :), rr_deflected(2, :), 'b-+', 'LineWidth', 2)
xlim([0. 0.001])
axis equal
grid on

%% functions --------------------------------------------------------------
function yq = hermite_spline(x, y, xq, tang_start, tang_end)

    % Input data
    % x:  x database
    % y:  y database
    % xq: x query
    % tang_start: tangent coefficient at x(1)
    % tang_end:   tangent coefficient at x(end)
    %
    % Output data
    % yq: y interpolated

    % build tangent vector using centered differences
    m = zeros(1, numel(x));
    m(1) = tang_start;
    m(end) = tang_end;
    % for i = 2:numel(x) - 1
    %     m(i) = (y(i - 1) - y(i))/(x(i - 1) - x(i));
    % end
    for i = 2:numel(x) - 1
        m(i) = 0.5 * ((y(i + 1) - y(i)) / (x(i + 1) - x(i)) + ...
            (y(i) - y(i - 1)) / (x(i) - x(i - 1)));
    end

    yq = zeros(1, numel(xq));

    for i = 1:numel(xq)

        for j = 1:numel(x) - 1
            % two consecutive points x for which the point xq belongs
            if (xq(i) >= x(j)) && (xq(i) <= x(j + 1))
                t = (xq(i) - x(j)) / (x(j + 1) - x(j));

                h00 = hermite_p1(t);
                h10 = hermite_p2(t);
                h01 = hermite_p3(t);
                h11 = hermite_p4(t);

                yq(i) = h00 * y(j) + ...
                    h10 * (x(j + 1) - x(j)) * m(j) + ...
                    h01 * y(j + 1) + ...
                    h11 * (x(j + 1) - x(j)) * m(j + 1);
            end

        end

    end

end

function h = hermite_p1(t)
    h = 2 * t ^ 3 - 3 * t ^ 2 + 1;
end

function h = hermite_p2(t)
    h = t ^ 3 - 2 * t ^ 2 + t;
end

function h = hermite_p3(t)
    h = -2 * t ^ 3 + 3 * t ^ 2;
end

function h = hermite_p4(t)
    h = t ^ 3 - t ^ 2;
end

function [dcsi_te, dcsi_le] = geoseries(start_x, end_x, n_elem, r)

    dl = zeros(1, n_elem);
    dcsi_te = zeros(1, n_elem + 1);
    % check series convergence
    if r >= 1 || r <= 0
        error('insert a growth value between 0 and 1')
    end

    r = 1 - r; % growth ratio (intended as decreasing ratio)
    % get size of the first element knowing the sum of the geometric series
    dl(1) = (1 - r) / (1 - r ^ (n_elem));

    for i = 2:n_elem
        dl(i) = r * dl(i - 1); % geometric series-> get length of element
    end

    for i = 1:n_elem
        dcsi_te(i + 1) = dcsi_te(i) + dl(i); %-> position
    end

    dcsi_le = fliplr(1 - dcsi_te);
    % rescale in the interval
    dcsi_te = rescale_csi(dcsi_te, start_x, end_x);
    dcsi_le = rescale_csi(dcsi_le, start_x, end_x);
end

function dcsi = geoseries_both(start_x, end_x, n_elem, r_le, r_te)

    n_elem_le = ceil(n_elem / 2);
    n_elem_te = floor(n_elem / 2);
    mid_point = (start_x + end_x) / 2;
    start_x_le = start_x;
    end_x_le = mid_point;
    start_x_te = mid_point;
    end_x_te = end_x;
    [~, dcsi_le] = geoseries(start_x_le, end_x_le, n_elem_le, r_le);
    [dcsi_te, ~] = geoseries(start_x_te, end_x_te, n_elem_te, r_te);
    dcsi = [dcsi_le(1:end - 1) dcsi_te];
end

function dcsi = geoseries_hinge(start_x, end_x, n_elem, xh, ...
        r_le_aft, r_te_aft, r_le_fore, r_te_fore)

    n_elem_aft = ceil(n_elem * xh);
    n_elem_fore = n_elem - n_elem_aft;
    start_x_aft = start_x;
    end_x_aft = xh;
    start_x_fore = xh;
    end_x_fore = end_x;
    dcsi_aft = geoseries_both(start_x_aft, end_x_aft, n_elem_aft, r_le_aft, r_te_aft);
    dcsi_fore = geoseries_both(start_x_fore, end_x_fore, n_elem_fore, r_le_fore, r_te_fore);
    dcsi = [dcsi_aft dcsi_fore];

end

function dcsi_rescaled = rescale_csi(dcsi, rmin, rmax)
    dcsi_rescaled = dcsi * (rmax - rmin) + rmin;
end

% test hinge deflection
function rr = hinge_deflection(rr, theta, xh, offset)
    position = [xh; 0];

    for ip = 1:size(rr, 2)
        d = (rr(:, ip) - position)' * [1; 0];

        if d >= offset
            Rot = [cos(theta), -sin(theta); ...
                       sin(theta), cos(theta)];
            rr(:, ip) = position + Rot * (rr(:, ip) - position);
        elseif d >= -offset
            u = offset;
            m = -cot(theta);
            yc = -m * u * (1 + cos(theta)) + u * sin(theta);
            thp = 0.5 * (d + u) / u * theta;
            rr(:, ip) = [yc * sin(thp) - u + position(1) - rr(2, ip) * sin(thp); ...
                            yc * (1 .- cos(thp)) + rr(2, ip) * cos(thp)];
        end

    end

end

%% %%-- just for verification --%%
function [x, y] = naca4digit(M, P, SS, c, n)

    m = M / 100;
    p = P / 10;
    t = SS / 100;

    if (m == 0) % p must be .ne. 0
        p = 1;
    end

    %> === Chord discretization ===
    xv = linspace(0.0, c, n + 1); % uniform
    % xv = c/2.0 .*(1.0-cos(pi.*xv./c));   % cosine
    xv = c .* (1.0 - cos(0.5 * pi .* xv ./ c)); % half-cosine

    %> === Thickness ===
    ytfcn = @(x) 5 .* t .* c .* (0.2969 .* (x ./ c) .^ 0.5 - 0.1260 .* (x ./ c) ...
        - 0.3516 .* (x ./ c) .^ 2 + 0.2843 .* (x ./ c) .^ 3 - 0.1015 .* (x ./ c) .^ 4);
    % 1015: open TE
    % 1036: closed TE
    yt = ytfcn(xv);

    %> === Mean line ===
    yc = zeros(size(xv));

    for ii = 1:n + 1

        if xv(ii) <= p * c
            yc(ii) = c * (m / p ^ 2 * (xv(ii) / c) * (2 * p - (xv(ii) / c)));
        else
            yc(ii) = c * (m / (1 - p) ^ 2 * (1 + (2 * p - (xv(ii) / c)) * (xv(ii) / c) -2 * p));
        end

    end

    %> === Mean line slope ===
    dyc = zeros(size(xv));

    for ii = 1:n + 1

        if xv(ii) <= p * c
            dyc(ii) = m / p ^ 2 * 2 * (p - xv(ii) / c);
        else
            dyc(ii) = m / (1 - p) ^ 2 * 2 * (p - xv(ii) / c);
        end

    end

    %> === Upper and lower sides ===
    th = atan2(dyc, 1);
    xU = xv - yt .* sin(th);
    yU = yc + yt .* cos(th);
    xL = xv + yt .* sin(th);
    yL = yc - yt .* cos(th);

    %> Sort points
    x = zeros(1, 2 * n + 1);
    y = zeros(1, 2 * n + 1);

    for ii = 1:n
        x(ii) = xL(n + 2 - ii);
        y(ii) = yL(n + 2 - ii);
    end

    x(n + 1:2 * n + 1) = xU;
    y(n + 1:2 * n + 1) = yU;
end

function chain_points = catmull_rom_chain(points, num_points)
    % Calculate Catmull-Rom for a sequence of initial points and return the combined curve.
    % :param points: Base points from which the quadruples for the algorithm are taken
    % :param num_points: The number of points to include in each curve segment
    % :return: The chain of all points (points of all segments)

    QUADRUPLE_SIZE = 4;
    point_quadruples = zeros(4, 2, num_segments(points));

    for idx_segment_start = 1:num_segments(points)
        point_quadruples(:, :, idx_segment_start) = points(idx_segment_start:idx_segment_start + QUADRUPLE_SIZE - 1, :);
    end

    all_splines = zeros(num_points * num_segments(points), 2);
    j = 1;

    for i = 1:num_segments(points)
        all_splines(j:j + num_points - 1, :) = catmull_rom_spline(point_quadruples(1, :, i), point_quadruples(2, :, i), point_quadruples(3, :, i), point_quadruples(4, :, i), num_points);
        j = j + num_points;
    end

    chain_points = all_splines;
end

function num_segments = num_segments(points)
    num_segments = size(points, 1) - 3;
end

function spline_points = catmull_rom_spline(P0, P1, P2, P3, num_points)
    % Compute the points in the spline segment
    % :param P0, P1, P2, and P3: The (x,y) point pairs that define the Catmull-Rom spline
    % :param num_points: The number of points to include in the resulting curve segment
    % :param alpha: 0.5 for the centripetal spline, 0.0 for the uniform spline, 1.0 for the chordal spline.
    % :return: The points

    % Calculate t0 to t4. Then only calculate points between P1 and P2.
    % Reshape linspace so that we can multiply by the points P0 to P3
    % and get a point for each value of t.
    function t = tj(ti, pi, pj)
        alpha = 0.5;
        xi = pi(1);
        yi = pi(2);
        xj = pj(1);
        yj = pj(2);
        dx = xj - xi;
        dy = yj - yi;
        l = sqrt(dx ^ 2 + dy ^ 2);
        t = ti + l ^ alpha;
    end

    t0 = 0.0;
    t1 = tj(t0, P0, P1);
    t2 = tj(t1, P1, P2);
    t3 = tj(t2, P2, P3);
    t = linspace(t1, t2, num_points)';

    A1 = (t1 - t) / (t1 - t0) * P0 + (t - t0) / (t1 - t0) * P1;
    A2 = (t2 - t) / (t2 - t1) * P1 + (t - t1) / (t2 - t1) * P2;
    A3 = (t3 - t) / (t3 - t2) * P2 + (t - t2) / (t3 - t2) * P3;
    B1 = ((t2 - t) / (t2 - t0)) .* A1 + ((t - t0) / (t2 - t0)) .* A2;
    B2 = ((t3 - t) / (t3 - t1)) .* A2 + ((t - t1) / (t3 - t1)) .* A3;
    spline_points = [(t2 - t) / (t2 - t1) .* B1(:, 1) + (t - t1) / (t2 - t1) .* B2(:, 1), (t2 - t) / (t2 - t1) .* B1(:, 2) + (t - t1) / (t2 - t1) .* B2(:, 2)];
end
