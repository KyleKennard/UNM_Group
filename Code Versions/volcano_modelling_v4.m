clear all; close all; format long;

%% Import data              

volcano = 'Mayon';

dted_file = dir(strcat(volcano, '/*.dt2'));
[elevs, refvec] = dted(fullfile(dted_file.folder, dted_file.name));

latlim = [refvec(2)-1, refvec(2)];
lonlim = [refvec(3), refvec(3)+1];

% From Grosse database (lat/lon @centre, basal radius, height)

dat_file = dir(strcat(volcano, '/*.dat'));
volc_dat = importdata(fullfile(dat_file.folder, dat_file.name));

volc_lat = volc_dat(1);
volc_lon = volc_dat(2);
base_rad = volc_dat(6)/2;
height   = volc_dat(11);

line_adj = 1.25;   % frac. change from Rbase for extending profile lines

%% Define lat/lon lines over which to gather elevation data

ddeg = km2deg(base_rad*line_adj);

theta = (0:10:180);

px_arr = zeros(length(theta), 2);
py_arr = zeros(length(theta), 2);
pz_arr = zeros(length(theta), 2);

xline_arr = cell(1, length(theta));
yline_arr = cell(1, length(theta));

zline_arr = cell(1, length(theta));
range_arr = cell(1, length(theta));

for t = 1:length(theta)
    
    plat = volc_lat + ddeg*[sind(theta(t)), -sind(theta(t))];
    plon = volc_lon + ddeg*[cosd(theta(t)), -cosd(theta(t))];

    [zline, range, latline, lonline] = mapprofile(elevs, refvec, plat, plon);
    
    px_arr(t, :) = deg2km(plon);
    py_arr(t, :) = deg2km(plat);
    pz_arr(t, :) = [zline(1), zline(end)];
    
    xline_arr{t} = deg2km(lonline);
    yline_arr{t} = deg2km(latline);
    
    zline_arr{t} = zline;
    range_arr{t} = deg2km(range - mean(range));
    
%     figure(1); hold on;
%     worldmap(latlim, lonlim)
%     meshm(elevs, refvec, size(elevs))
%     plot3m(latline, lonline, zline, 'w', 'LineWidth', 2)
%     demcmap(elevs)

end

%% Testing linear trend subtraction method

[n, ~, p] = affine_fit([px_arr(:), py_arr(:), pz_arr(:)]);

[px, py] = meshgrid(px_arr(:), py_arr(:));

zpln = -(n(1)/n(3)*px + n(2)/n(3)*py - dot(n, p)/n(3));
zfit = -(n(1)/n(3)*px_arr(:) + n(2)/n(3)*py_arr(:) - dot(n, p)/n(3));
zadj = pz_arr(:) - zfit;

figure(2); hold on;
set(gca, 'FontSize', 18)
surf(px*(1e-3), py*(1e-3), zpln, 'facecolor', 'y', 'facealpha', 0.5)
surf(px*(1e-3), py*(1e-3), zeros(size(px)), 'facecolor', 'g', 'facealpha', 0.5)
plot3(px_arr(:)*(1e-3), py_arr(:)*(1e-3), pz_arr(:), 'ko', 'markerfacecolor', 'r');
plot3(px_arr(:)*(1e-3), py_arr(:)*(1e-3), zadj, 'ko', 'markerfacecolor', 'b');
xlabel('X (km)')
ylabel('Y (km)')
zlabel('Elevation (m)')
view(25, 5)
grid on;

%% Subtract any linear trend from elevations

for t = 1:length(theta)
    
    zt = zline_arr{t};
    rt = range_arr{t};
    
    zfit = -(n(1)/n(3)*xline_arr{t} + n(2)/n(3)*yline_arr{t} - dot(n, p)/n(3));
    zadj = zt - zfit;
    
    posz = (zadj > 0);
    zadj = zadj(posz);
    zt   = zt(posz);
    rt   = rt(posz);
    
    zline_arr{t} = zadj;
    range_arr{t} = rt;
    
    figure(3); hold on;
    set(gca, 'FontSize', 18)
    plot(rt, zt)
    plot(rt, zadj)
    xlabel('Range (km)')
    ylabel('Elevation (m)')
    grid on
    
end

%% Calculate and plot average profile elevations

rbins = 10^floor(log10(max(cellfun(@length, range_arr))));

range_min = min(cellfun(@min, range_arr));
range_max = max(cellfun(@max, range_arr));
range_bin = linspace(range_min, range_max, rbins);

zline_avg = zeros(1, length(range_bin)-1);
zline_err = zeros(1, length(range_bin)-1);
range_avg = zeros(1, length(range_bin)-1);
  
for r = 2:rbins
    
   zline_store = [];
   range_store = [];
  
   for t = 1:length(theta) 
       
       zl = zline_arr{t};
       rl = range_arr{t};
       
       Rinds = find(rl <= range_bin(r) & rl > range_bin(r-1));
       
       zline_store = [zline_store; zl(Rinds)];
       range_store = [range_store; rl(Rinds)];
       
   end
   
   zline_avg(r-1) = mean(zline_store); 
   zline_err(r-1) = std(zline_store);
   range_avg(r-1) = mean(range_store);
   
end

zline_avg = zline_avg(isfinite(zline_avg));
zline_err = zline_err(isfinite(zline_err));
range_avg = range_avg(isfinite(range_avg));

figure(3); hold on;
errorbar(range_avg, zline_avg, zline_err, 'k');

%% Build model from data

Rbase   = max(range_avg)*(1e3);
hmax    = max(zline_avg) - min(zline_avg);
tan_phi = tand(30);

r = linspace(0, Rbase, length(range_avg));

rcrit = (Rbase^2/hmax)*(tan_phi - sqrt(tan_phi^2 - (hmax/Rbase)^2));
alpha = (hmax - rcrit*tan_phi)*(rcrit*Rbase^2/(Rbase^2 - rcrit^2));

hr = zeros(size(r));

for i = 1:length(r)
    if r(i) < rcrit
        hr(i) = (rcrit - r(i))*tan_phi + alpha*(1/rcrit - rcrit/Rbase^2);
    else
        hr(i) = alpha*(1/r(i) - r(i)/Rbase^2);
    end
end

hr_plot  = [fliplr(hr) hr];
distance = linspace(-Rbase, Rbase, length(hr_plot))*(1e-3);

%% Plotting model with data 

figure(4); hold on;
set(gca, 'FontSize', 18)
plot(distance, hr_plot, 'Linewidth', 2)
errorbar(range_avg, zline_avg - min(zline_avg), zline_err);
xlabel('Distance (km)')
ylabel('Edifice Height (m)')
legend('model', 'data (with error)')
ylim([-100, max(zline_avg)])
grid on

%% Goodness-of-fit

%Calculate RMS between model and data

idx=knnsearch(distance',range_avg');
model_avg = hr_plot(idx);

RMS = rms(model_avg-(zline_avg-min(zline_avg)));

fprintf('RMS: %.3f \n',RMS)
