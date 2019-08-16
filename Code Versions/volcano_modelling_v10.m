clear all; close all; format long;

%% User input

region  = 'Halmahera-Indonesia';
volcano = 'Gamkonora';

line_adj = 1.0;   % frac. change from Rbase for extending profile lines

%% Determine .dted file and .dat file paths for given volcano

region_dir = dir(region);

in_reg_dir = {region_dir.name};                             % Pick out directories inside region
in_reg_dir = in_reg_dir([region_dir.isdir]);                % directory (omitting '.' and '..')
in_reg_dir = in_reg_dir(cellfun(@length, in_reg_dir) > 2);

if length(dir(fullfile(region, '*.dt2'))) == 1
    
    dted_file = dir(fullfile(region, '*.dt2'));
    
    for d = 1:length(in_reg_dir)
       
        volctype_path = fullfile(region, in_reg_dir{d});
        volctype_dir  = dir(volctype_path);
        
        if ~isempty(nonzeros(cellfun(@(c) strcmp(c, volcano), {volctype_dir.name})))           
            dat_file = dir(fullfile(volctype_path, volcano, strcat(volcano, '.dat'))); 
        end       
    end
    
else
    for d = 1:length(in_reg_dir)
       
        subreg_dir_path = fullfile(region, in_reg_dir{d});
        subreg_dir      = dir(subreg_dir_path);
        
        in_subreg_dir   = {subreg_dir.name};
        in_subreg_dir   = in_subreg_dir([subreg_dir.isdir]);
        in_subreg_dir   = in_subreg_dir(cellfun(@length, in_subreg_dir) > 2);
        
        for sd = 1:length(in_subreg_dir)
        
            volctype_path = fullfile(subreg_dir_path, in_subreg_dir{sd});
            volctype_dir  = dir(volctype_path);
        
            if ~isempty(nonzeros(cellfun(@(c) strcmp(c, volcano), {volctype_dir.name})))           
                dat_file  = dir(fullfile(volctype_path, volcano, strcat(volcano, '.dat'))); 
                dted_file = dir(fullfile(subreg_dir_path, '*.dt2'));
            end 
        end
    end
end

%% Import data
%  (1) elevs and lat/lon limits from dted file 
%  (2) from Grosse database: lat/lon at volcano centre, basal radius, height

[elevs, refvec] = dted(fullfile(dted_file.folder, dted_file.name));

latlim = [refvec(2)-1, refvec(2)];
lonlim = [refvec(3), refvec(3)+1];

volc_dat = importdata(fullfile(dat_file.folder, dat_file.name));

volc_lat = volc_dat(1);
volc_lon = volc_dat(2);
base_rad = volc_dat(6)/2;
height   = volc_dat(11);

%% Define lat/lon lines over which to gather elevation data

ddeg = km2deg(base_rad*line_adj);

theta = (0:1:180);

px_arr = zeros(length(theta), 2);
py_arr = zeros(length(theta), 2);
pz_arr = zeros(length(theta), 2);

xline_arr = cell(1, length(theta));
yline_arr = cell(1, length(theta));

zline_arr = cell(1, length(theta));
range_arr = cell(1, length(theta));

% figure(1); hold on;
% worldmap(latlim, lonlim)
% meshm(elevs, refvec, size(elevs))
% demcmap(elevs)

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
%     plot3m(latline, lonline, zline, 'w')

end

%% Testing linear trend subtraction method

pXYZ = [px_arr(:), py_arr(:), pz_arr(:)];

[n, p]   = bestfit_plane(pXYZ);
[px, py] = meshgrid(pXYZ(:,1), pXYZ(:,2));

zpln = (dot(n, p) - n(1)*px - n(2)*py)/n(3);
zfit = (dot(n, p) - n(1)*pXYZ(:,1) - n(2)*pXYZ(:,2))/n(3);
zadj = pXYZ(:,3) - zfit;

figure(2); hold on;
set(gca, 'FontSize', 18)
surf(px*(1e-3), py*(1e-3), zpln*(1e-3), 'facecolor', 'y', 'facealpha', 0.5)
surf(px*(1e-3), py*(1e-3), zeros(size(px)), 'facecolor', 'g', 'facealpha', 0.5)
plot3(pXYZ(:,1)*(1e-3), pXYZ(:,2)*(1e-3), pXYZ(:,3)*(1e-3), 'ko', 'markerfacecolor', 'r');
plot3(pXYZ(:,1)*(1e-3), pXYZ(:,2)*(1e-3), zadj*(1e-3), 'ko', 'markerfacecolor', 'b');
xlabel('X (km)')
ylabel('Y (km)')
zlabel('Elevation (km)')
view(25, 5)
grid on;

%% Subtract any linear trend from elevations

for t = 1:length(theta)
    
    zt = zline_arr{t};
    rt = range_arr{t};

    zfit = (dot(n, p) - n(1)*xline_arr{t} - n(2)*yline_arr{t})/n(3);
    zadj = zt - zfit;
    
    posz = (zadj > 0);
    zadj = zadj(posz);
    zt   = zt(posz);
    rt   = rt(posz);
    
    zline_arr{t} = zadj - min(zadj);
    range_arr{t} = rt;
    
    figure(3); hold on;
    set(gca, 'FontSize', 18)
%     plot(rt, zt)
    plot(rt, zadj)
    xlabel('Range (km)')
    ylabel('Elevation (m)')
    grid on
    
end

%% Calculate and plot average profile elevations

window = 0.3;       % Optimized for Mayon: window = 0.3, step = 0.041
step   = 0.041;

range_min = min(cellfun(@min, range_arr)) - window;
range_max = max(cellfun(@max, range_arr));
range_bin = (range_min:step:range_max);
left      = range_bin(1);

high_flag = [];
cut_off   = 1.5;      % Height at which to determine outlier profiles

for r = 1:length(range_bin)
    
   right = left + window;
    
   zline_store = [];
   range_store = [];
  
   for t = 1:length(theta) 
       
       zl    = zline_arr{t};
       rl    = range_arr{t};      
       Rmask = (rl < right & rl >= left);
       
       zline_store = [zline_store; zl(Rmask)];
       range_store = [range_store; rl(Rmask)];
       
   end
   
   avgstore = mean(zline_store);
   
   for t = 1:length(theta)
       
       zl    = zline_arr{t};
       rl    = range_arr{t};      
       Rmask = (rl < right & rl >= left);

       q25 = prctile(zline_store, 25);
       q75 = prctile(zline_store, 75);
       
       iqr = q75 - q25;
       iqr_cutoff = iqr*cut_off;
       upper = avgstore + iqr_cutoff;
       
       check = zl(Rmask) > upper;
       
       if ~isempty(nonzeros(check))
           high_flag = [high_flag; t];
       end
   end 
   
   left = left + step;
   
end

%% Remove "too anomalous" profiles from average and errors

tarr  = (1:1:t);
tmask = ~ismember(tarr, tarr(unique(high_flag)));

zline_arr_cut = zline_arr(tmask);
range_arr_cut = range_arr(tmask);

[zline_avg, zline_err, range_avg] = deal([]);

left = range_bin(1);

for r = 1:length(range_bin)
    
   right = left + window;
    
   zline_store = [];
   range_store = [];
   
   for t = 1:length(zline_arr_cut)
       
       zl = zline_arr_cut{t};
       rl = range_arr_cut{t};
       
       if r == 1
           figure(4); hold on;
           set(gca, 'FontSize', 18)
           plot(rl, zl)
           xlabel('Range (km)')
           ylabel('Elevation (m)')
           grid on
       end
       
       Rmask = (rl < right & rl >= left);
       
       zline_store = [zline_store; zl(Rmask)];
       range_store = [range_store; rl(Rmask)];
   end
   
   if length(zline_store) >= 2
       
       zline_avg = [zline_avg, mean(zline_store)];
       zline_err = [zline_err, std(zline_store)];
       range_avg = [range_avg, mean(range_store)];       
   end  
   
   left = left + step;
   
end

figure(4); hold on; plot(range_avg, zline_avg, 'k', 'linewidth', 2);

%% Build model from data

Rbase   = max(range_avg)*(1e3);
hmax    = max(zline_avg);
tan_phi = tand(30);

rmod  = linspace(0, Rbase, length(range_avg));
rcrit = (Rbase^2/hmax)*(tan_phi - sqrt(tan_phi^2 - (hmax/Rbase)^2));
alpha = (hmax - rcrit*tan_phi)*(rcrit*Rbase^2/(Rbase^2 - rcrit^2));

hr = zeros(size(rmod));

for i = 1:length(rmod)
    if rmod(i) < rcrit
        hr(i) = (rcrit - rmod(i))*tan_phi + alpha*(1/rcrit - rcrit/Rbase^2);
    else
        hr(i) = alpha*(1/rmod(i) - rmod(i)/Rbase^2);
    end
end

hr_plot  = [fliplr(hr) hr];
distance = linspace(-Rbase, Rbase, length(hr_plot))*(1e-3);

%% Goodness of Fit (calculates RMS between model and data)

idx       = knnsearch(distance', range_avg');
model_avg = hr_plot(idx);

residual  = zline_avg - model_avg;
chi_sq    = sum((residual./zline_err).^2)/(length(residual)-1);

fprintf('chi^2: %.5f \n', chi_sq)

%% Plotting model with data 

figtitle = sprintf('%s || \\chi^2_{\\nu} = %.5f || lineadj = %.1f', volcano, chi_sq, line_adj);

figure(5); hold on;
set(gca, 'FontSize', 16)

plot(distance, hr_plot, 'Linewidth', 2)
errorbar(range_avg, zline_avg, zline_err);
xlabel('Distance (km)')
ylabel('Edifice Height (m)')
title(figtitle,'interpreter','tex')
legend('model', 'data (with error)')
set(legend, 'FontSize', 16)
ylim([-10, max(zline_avg) + 100])
grid on

%% Save volcano-model fit figure

% sv_fig_file = sprintf('%s_chisq%.5f_lineadj%.1f.fig', volcano, chi_sq, line_adj);
% sv_fig_path = fullfile(dat_file.folder, sv_fig_file);
% savefig(5, sv_fig_path);
% 
% %% Save goodness of fit in designated .dat file
% 
% sv_dat_file = strcat(volcano, '_goodfit.dat');
% sv_dat_path = fullfile(dat_file.folder, sv_dat_file); 
% 
% if isfile(sv_dat_path)   
%     good_fit = importdata(sv_dat_path);
%     good_fit(:, end+1) = [line_adj, chi_sq];
%     save(sv_dat_path, 'good_fit', '-ascii');   
% else    
%     good_fit = [line_adj; chi_sq];
%     save(sv_dat_path, 'good_fit', '-ascii');
% end

%% Plane of best fit function. Input is an M x 3 array of x,y,z points.
%  Output is a normal vector to the corresponding plane of best fit, and a
%  point in that plane.

function [normal, point] = bestfit_plane(XYZarr)

    point   = mean(XYZarr);
    redmat  = XYZarr - point;
    [~,~,V] = svd(redmat);
    
    if acosd(V(3,1)) > 0.5 && acosd(V(3,1)) < 178.0   
        normal = V(:,3);
    else                     
        normal = V(:,1);
    end
  
end
