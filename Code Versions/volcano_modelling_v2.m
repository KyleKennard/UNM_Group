clear all; close all; format long;

%% Import data              

volcano = 'Mayon';

dted_file = dir(strcat(volcano, '/*.dt2'));
[elevs, refvec] = dted(fullfile('Mayon', 'n13_e123_1arc_v3.dt2'));

latlim = [refvec(2)-1, refvec(2)];
lonlim = [refvec(3), refvec(3)+1];

% From Grosse database (lat/lon @centre, basal radius, height)

dat_file = dir(strcat(volcano, '/*.dat'));
volc_dat = importdata(fullfile('Mayon', 'Mayon.dat'));

volc_lat = volc_dat(1);
volc_lon = volc_dat(2);
base_rad = volc_dat(6)/2;
height   = volc_dat(11);

line_adj = 1.5;   % frac. change from Rbase for extending profile lines

%% Define lat/lon lines over which to gather elevation data

ddeg = km2deg(base_rad*line_adj);

theta = (0:15:180);

latline_arr = cell(1, length(theta));
lonline_arr = cell(1, length(theta));
zline_arr   = cell(1, length(theta));
range_arr   = cell(1, length(theta));

for t = 1:length(theta)
    
    plat = volc_lat + ddeg*[sind(theta(t)), -sind(theta(t))];
    plon = volc_lon + ddeg*[cosd(theta(t)), -cosd(theta(t))];

    [zline, range, latline, lonline] = mapprofile(elevs, refvec, plat, plon);
    
    latline_arr{t} = latline;
    lonline_arr{t} = lonline;
    zline_arr{t}   = zline;
    range_arr{t}   = deg2km(range - mean(range));
    
    % Plotting
    
%     figure(1); hold on;
%     worldmap(latlim, lonlim)
%     meshm(elevs, refvec, size(elevs))
%     plot3m(latline, lonline, zline, 'w', 'LineWidth', 2)
%     demcmap(elevs)
    
    figure(2); hold on;
    set(gca, 'FontSize', 18)
    plot(deg2km(range - mean(range)), zline)
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
range_avg = zeros(1, length(range_bin)-1);

zlinestores=cell(rbins-1,1);
  
for r = 2:rbins
    
   nterms    = 0;
   zline_sum = 0;
   range_sum = 0;
  zline_storet=zeros(max(cellfun(@length, range_arr)),max(cellfun(@length, range_arr)));
   for t = 1:length(theta) 
       
       zl = zline_arr{t};
       rl = range_arr{t};
       
       Rinds = find(rl <= range_bin(r) & rl > range_bin(r-1));
       
       zline_sum = zline_sum + sum(zl(Rinds));
       range_sum = range_sum + sum(rl(Rinds));
       nterms    = nterms + length(Rinds); 
       
       zline_storet(t,1:length(Rinds)) = zl(Rinds);
   end
   
   zline_storet=zline_storet(:);
   zline_storet=zline_storet(find(zline_storet~=0));
   
   zline_storer{r-1} = zline_storet;
   
   zline_avg(r-1) = zline_sum/nterms; 
   range_avg(r-1) = range_sum/nterms;
end

zline_avg  = zline_avg(isfinite(zline_avg));
range_avg = range_avg(isfinite(range_avg));

figure(2); hold on;
plot(range_avg, zline_avg, 'k', 'LineWidth', 2)
%% Standard Deviation
for r=2:rbins
    err(r-1)=std(zline_storer{r-1});
end
%% Build model from data

Rbase   = max(range_avg)*(1e3);
hmax    = max(zline_avg) - min(zline_avg);
tan_phi = tand(30);

r = linspace(0, Rbase, 1e3);

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

figure(3); hold on;
set(gca, 'FontSize', 18)
plot(distance, hr_plot, '--')
plot(range_avg, zline_avg - min(zline_avg))
xlabel('Distance (km)')
ylabel('Edifice Height (m)')
legend('model', 'data')
grid on

%% Plotting Error Bars
figure(4);
set(gca, 'FontSize', 18)
errorbar(range_avg, zline_avg-min(zline_avg),err,'k');
xlabel('Range (km)')
ylabel('Elevation (m)')
grid on
