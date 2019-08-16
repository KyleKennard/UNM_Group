clear all; close all; format long;

%% User input

region  = 'Halmahera-Indonesia';
volcano = 'Tidore';

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

theta = (0:1:359);

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
    
    plat =[volc_lat, volc_lat + ddeg*sind(theta(t))];
    plon = [volc_lon, volc_lon + ddeg*cosd(theta(t))];

    [zline, range, latline, lonline] = mapprofile(elevs, refvec, plat, plon);
    
    px_arr(t, :) = deg2km(plon);
    py_arr(t, :) = deg2km(plat);
    pz_arr(t, :) = [zline(1), zline(end)];
    
    xline_arr{t} = deg2km(lonline);
    yline_arr{t} = deg2km(latline);
    
    zline_arr{t} = zline;
    if   t <= 90 | t > 270
        range_arr{t} = deg2km(range);
    elseif  t <= 270 & t > 90
        range_arr{t}=deg2km(-range);
    end
%     
%     figure(1); hold on;
%     plot3m(latline, lonline, zline, 'w')

end

%% Get user input as to which normal vector to use in planar fit correction
pXYZ = [px_arr(:,2), py_arr(:,2), pz_arr(:,2)];
[n1, n3, p] = bestfit_plane(pXYZ);
[px, py] = meshgrid(pXYZ(:,1), pXYZ(:,2));

% Plane fit using n1 = V(:,1) normal vector

zpln1 = (dot(n1, p) - n1(1)*px - n1(2)*py)/n1(3);
zfit1 = (dot(n1, p) - n1(1)*pXYZ(:,1) - n1(2)*pXYZ(:,2))/n1(3);
zadj1 = pXYZ(:,3) - zfit1;

% Plot a sample planar fit

figure(2); hold on;
set(gca, 'FontSize', 14)

hold on;
set(gca, 'FontSize', 14)
surf(px*(1e-3), py*(1e-3), zeros(size(px)), 'facecolor', 'g', 'facealpha', 0.8, 'edgecolor', 'none');
plot3(pXYZ(:,1)*(1e-3), pXYZ(:,2)*(1e-3), zadj1*(1e-3), 'ko', 'markerfacecolor', 'b');
xlabel('X (km)');   ylabel('Y (km)');   zlabel('Elevation (km)');
title('Sample Plane Fit')
view(25, 5)
grid on;

%Plot a circle that will show angles in increments of 10 degrees
tick=linspace(10,360,36);
for q=1:length(tick)
    tickline=[base_rad/1000*cosd(tick(q)),base_rad/1000*sind(tick(q)),0];
    plot3([xline_arr{1}(1)/1000,xline_arr{1}(1)/1000+tickline(1)],[yline_arr{1}(1)/1000,yline_arr{1}(1)/1000+tickline(2)],[0,tickline(3)],'k');
end
plot3([xline_arr{1}(1)/1000,xline_arr{1}(1)/1000+base_rad/1000],[yline_arr{1}(1)/1000,yline_arr{1}(1)/1000],[0,0],'m')
   
    
%Get user to decide which azimuths to exclude

excludes=input('Input a matrix of azimuth angles to exlude from processing. \n');
theta=setdiff(theta,excludes);
close all
%% Redo Planar Fit with chosen azimuths excluded
if ~isempty(excludes)
    pln1=[]; zfit1=[]; zadj1=[]; zpln3=[]; zfit3=[]; zadj3=[]; n1=[]; n3=[]; p=[]; px=[]; py=[];
    pXYZ=zeros(length(theta),3);
    for t=1:length(theta)
        pXYZ(t,:)=[px_arr(theta(t),2),py_arr(theta(t),2),pz_arr(theta(t),2)];
    end
end
    [n1, n3, p] = bestfit_plane(pXYZ);
    [px, py] = meshgrid(pXYZ(:,1), pXYZ(:,2));

% Plane fit using n1 = V(:,1) normal vector

zpln1 = (dot(n1, p) - n1(1)*px - n1(2)*py)/n1(3);
zfit1 = (dot(n1, p) - n1(1)*pXYZ(:,1) - n1(2)*pXYZ(:,2))/n1(3);
zadj1 = pXYZ(:,3) - zfit1;

% Plane fit using n3 = V(:,3) normal vector

zpln3 = (dot(n3, p) - n3(1)*px - n3(2)*py)/n3(3);
zfit3 = (dot(n3, p) - n3(1)*pXYZ(:,1) - n3(2)*pXYZ(:,2))/n3(3);
zadj3 = pXYZ(:,3) - zfit3;

% Plotting

figure(2); hold on;
set(gca, 'FontSize', 14)

subplot(2,1,1); hold on;
set(gca, 'FontSize', 14)
surf(px*(1e-3), py*(1e-3), zpln1*(1e-3), 'facecolor', 'y', 'facealpha', 0.8, 'edgecolor', 'none');
surf(px*(1e-3), py*(1e-3), zeros(size(px)), 'facecolor', 'g', 'facealpha', 0.8, 'edgecolor', 'none');
plot3(pXYZ(:,1)*(1e-3), pXYZ(:,2)*(1e-3), pXYZ(:,3)*(1e-3), 'ko', 'markerfacecolor', 'r');
plot3(pXYZ(:,1)*(1e-3), pXYZ(:,2)*(1e-3), zadj1*(1e-3), 'ko', 'markerfacecolor', 'b');
xlabel('X (km)');   ylabel('Y (km)');   zlabel('Elevation (km)');
title('normal vector: n1')
view(25, 5)
grid on;

subplot(2,1,2); hold on;
set(gca, 'FontSize', 14)
surf(px*(1e-3), py*(1e-3), zpln3*(1e-3), 'facecolor', 'y', 'facealpha', 0.8, 'edgecolor', 'none');
surf(px*(1e-3), py*(1e-3), zeros(size(px)), 'facecolor', 'g', 'facealpha', 0.8, 'edgecolor', 'none');
plot3(pXYZ(:,1)*(1e-3), pXYZ(:,2)*(1e-3), pXYZ(:,3)*(1e-3), 'ko', 'markerfacecolor', 'r');
plot3(pXYZ(:,1)*(1e-3), pXYZ(:,2)*(1e-3), zadj3*(1e-3), 'ko', 'markerfacecolor', 'b');
xlabel('X (km)');   ylabel('Y (km)');   zlabel('Elevation (km)');
title('normal vector: n3')
view(25, 5)
grid on;

%User Input
n13 = input('Which normal vector would you like to use? Please enter either 1 or 3.\n');

if      n13 == 1;    n = n1;
elseif  n13 == 3;    n = n3;
else;   error('Must enter either 1 or 3 at input prompt.');
end

%% Subtract any linear trend from elevations

for t = 1:length(theta)
    
        zt = zline_arr{theta(t)};
        rt = range_arr{theta(t)};

        zfit = (dot(n, p) - n(1)*xline_arr{theta(t)} - n(2)*yline_arr{theta(t)})/n(3);
        zadj = zt - zfit;

        posz = (zadj > 0);
        zadj = zadj(posz);
        zt   = zt(posz);
        rt   = rt(posz);

        zline_arr{theta(t)} = zadj - min(zadj);
        range_arr{theta(t)} = rt;

        figure(3); hold on;
        set(gca, 'FontSize', 18)
        plot(rt, zadj)
        xlabel('Range (km)')
        ylabel('Elevation (m)')
        grid on
end

%% Calculate and plot average profile elevations

window = 11/30;       % Optimized for Mayon: window = 11/30, step = 4/75
step   = 4/75;

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

           zl    = zline_arr{theta(t)};
           rl    = range_arr{theta(t)};      
           Rmask = (rl < right & rl >= left);

           zline_store = [zline_store; zl(Rmask)];
           range_store = [range_store; rl(Rmask)];

   end
   
   avgstore = mean(zline_store);
   
   for t = 1:length(theta)

           zl    = zline_arr{theta(t)};
           rl    = range_arr{theta(t)};      
           Rmask = (rl < right & rl >= left);

           q25 = prctile(zline_store, 25);
           q75 = prctile(zline_store, 75);

           iqr = q75 - q25;
           iqr_cutoff = iqr*cut_off;
           upper = avgstore + iqr_cutoff;

           check = zl(Rmask) > upper;

           if ~isempty(nonzeros(check))
               high_flag = [high_flag; theta(t)];
           end

   end 
   
   left = left + step;
   
end

%% Remove "too anomalous" profiles from average and errors

tarr  = theta;
tmask = ~ismember(tarr, unique(high_flag));

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
       
   zline_avg = [zline_avg, mean(zline_store)];
   zline_err = [zline_err, std(zline_store)];
   range_avg = [range_avg, mean(range_store)];
  
   left = left + step;
  
end

figure(4); hold on; plot(range_avg, zline_avg, 'k', 'linewidth', 2);

%% Cut out zero or infinite or NaN values from arrays, if applicable

clean_err = (zline_err ~= 0 & isfinite(zline_err));
zline_avg = zline_avg(clean_err);
zline_err = zline_err(clean_err);
range_avg = range_avg(clean_err);

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

function [norm1, norm3, point] = bestfit_plane(XYZarr)

    point   = mean(XYZarr);
    redmat  = XYZarr - point;
    [~,~,V] = svd(redmat);
    
    norm1 = V(:,1);
    norm3 = V(:,3);
    
end
