clear all; close all; format long;

%% Import data              

dted_file = 'n13_e123_1arc_v3.dt2';

% From Grosse database (lat/lon of volcano centre)

volc_lat = 13.255;
volc_lon = 123.686;

basal_width= 8.89;
%% Define profile lines over which to get elevations

[elevs, refvec] = dted(dted_file);

latlim = [refvec(2)-1, refvec(2)];
lonlim = [refvec(3), refvec(3)+1];

ddeg = km2deg(basal_width/2);     % Determines how long the profile lines should be

volclat_min = volc_lat - ddeg;
volclat_max = volc_lat + ddeg;
volclon_min = volc_lon - ddeg;
volclon_max = volc_lon + ddeg;

plat1   = [volclat_min, volclat_max];
plon1   = [volclon_min, volclon_max];

plat2   = [volclat_max, volclat_min];
plon2   = [volclon_min, volclon_max];

plat3   = [volclat_min-(1e-4), volclat_max-(1e-4)];
plon3   = [mean(plon1), mean(plon1)];

plat4   = [mean(plat1), mean(plat1)];
plon4   = [volclon_min-(1e-4), volclon_max-(1e-4)];

[zline1, range1, latline1, lonline1] = mapprofile(elevs, refvec, plat1, plon1);
[zline2, range2, latline2, lonline2] = mapprofile(elevs, refvec, plat2, plon2);
[zline3, range3, latline3, lonline3] = mapprofile(elevs, refvec, plat3, plon3);
[zline4, range4, latline4, lonline4] = mapprofile(elevs, refvec, plat4, plon4);

zline_avg = (zline1 + zline2 + zline3 + zline4)/4;
range_avg = deg2km((range1 + range2 + range3 + range4)/4);

%% Error Bars

err=sqrt((1/3)*((zline1-zline_avg).^2+(zline2-zline_avg).^2+(zline3-zline_avg).^2+(zline4-zline_avg).^2));

%% Plotting

%figure(1); hold on;
%worldmap(latlim, lonlim)
%meshm(elevs, refvec, size(elevs))
%demcmap(elevs)
%plot3m(latline1, lonline1, zline1, 'w', 'LineWidth', 2)
%plot3m(latline2, lonline2, zline2, 'w', 'LineWidth', 2)
%plot3m(latline3, lonline3, zline3, 'w', 'LineWidth', 2)
%plot3m(latline4, lonline4, zline4, 'w', 'LineWidth', 2)

figure(2); hold on;
set(gca, 'FontSize', 18)
plot(deg2km(range1 - mean(range1)), zline1)
plot(deg2km(range2 - mean(range2)), zline2)
plot(deg2km(range3 - mean(range3)), zline3)
plot(deg2km(range4 - mean(range4)), zline4)
plot(range_avg - mean(range_avg), zline_avg, 'k', 'LineWidth', 2)
xlabel('Range (km)')
ylabel('Elevation (m)')
grid on

figure(3); hold on;
plot(deg2km(range1 - mean(range1)), zline1)
plot(deg2km(range2 - mean(range2)), zline2)
plot(deg2km(range3 - mean(range3)), zline3)
plot(deg2km(range4 - mean(range4)), zline4)
plot(range_avg - mean(range_avg), zline_avg, 'k', 'LineWidth', 2)
errorbar(range_avg(1:3:321) - mean(range_avg),zline_avg(1:3:321),err(1:3:321),'k')
xlabel('Range (km)')
ylabel('Elevation (m)')
grid on