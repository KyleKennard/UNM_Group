%% Calculate a given muon trajectory's path length through individual voxels

clear all; close all;

%% Define a detector location and look angle

det   = [1.0; 1.0; 1.0];    % x,y,z
theta = 50;                 % from zenith
phi   = 40;                 % asimuthal from north

%% Create the grid (distances are all in m)

x_max  = 10;
y_max  = x_max;
vside  = 0.5;

z_max  = 2;
z_min  = 0;
dz     = vside;

lambda = 2*sqrt(z_max)^2;         % for defining Gaussian hill

xarr   = 0:vside:x_max;
yarr   = 0:vside:y_max;

xshifted = wshift('1D', xarr, 1);
xshifted(end) = max(xarr) + vside;

yshifted = wshift('1D',yarr,1);
yshifted(end) = max(yarr) + vside;

[XI,YI] = meshgrid(xarr, yarr);
[XIshifted, YIshifted] = meshgrid(xshifted, yshifted);

% ElevI = max_z*ones(size(XI));
ElevI = 5*z_max*exp(-((XI - x_max/2).^2 + (YI - y_max/2).^2)/lambda) + z_max;

nlayers    = floor((ElevI-z_min)/dz);
maxlayers  = max(nlayers(:));
num_voxels = sum(nlayers(:));

for ilayer = 1:maxlayers
    
    layerbot = (z_min+dz*(ilayer-1))*ones(size(ElevI));
    layertop = (z_min+dz*ilayer)*ones(size(ElevI));
    allow_x  = XI(ElevI >= layertop);
    allow_x_shift = XIshifted(ElevI >= layertop);
    allow_y = YI(ElevI >= layertop);
    allow_y_shift = YIshifted(ElevI >= layertop);
    allow_bot  = layerbot(ElevI >= layertop);
    allow_top  = layertop(ElevI >= layertop);
    bot_corner = [allow_x(:)'; allow_y(:)'; allow_bot(:)'];
    top_corner = [allow_x_shift(:)'; allow_y_shift(:)'; allow_top(:)'];
    
    if ilayer == 1
        
        voxel_corner = bot_corner;
        voxel_top_corner = top_corner;
        
    else
        
        voxel_corner = [voxel_corner bot_corner];
        voxel_top_corner = [voxel_top_corner top_corner];
        
    end
end

zarr = z_min:dz:max(ElevI(:));

%% Draw trajectory and plot

% r = 10;
% 
% x = r*sind(theta)*sind(phi)+det(1);
% y = r*sind(theta)*cosd(phi)+det(2);
% z = r*cosd(theta)+det(3);
% 
% figure(); hold on;
% scatter3(voxel_corner(1,:),voxel_corner(2,:),voxel_corner(3,:),'b.');
% % scatter3(voxel_top_corner(1,:),voxel_top_corner(2,:),voxel_top_corner(3,:),'g.');
% plot3([det(1), x], [det(2), y], [det(3), z], 'r')
% plot3(det(1), det(2), det(3), 'kd', 'markerfacecolor', 'y', 'markersize', 10)
% grid on

%% Create a global index matrix

global_index = zeros(size(voxel_corner));

for i = 1:size(voxel_corner, 2)
    
    global_index(1,i) = find(abs(xarr - voxel_corner(1,i))==0);
    global_index(2,i) = find(abs(yarr - voxel_corner(2,i))==0);
    global_index(3,i) = find(abs(zarr - voxel_corner(3,i))==0);
    
end

%% Kat's codes

Lij_slice = Lij_create_v2(theta, phi, det, global_index, voxel_corner, xarr, yarr, zarr, ElevI, maxlayers, dz);

%% Remove 0 elements for Lij_slice

slice = Lij_slice(Lij_slice~=0);

avg = mean(slice);

figure()

bar(1:length(slice),slice)
hold on
plot(1:length(slice),ones(1,length(slice)).*avg,'r--')

title(['Ascension = ',num2str(phi),' deg, Declination = ',num2str(theta),' deg']);
xlabel('Voxel Number')
ylabel('Voxel Path Length (m)')

%% Draw voxels through which trajectory passes

% Define a colormap

figure()

for i=1:length(Lij_slice)
    if Lij_slice(i)~=0
        hits(:,i)=voxel_corner(:,i);
        Xbound=[hits(1,i);hits(1,i)+vside;hits(1,i)+vside;hits(1,i);hits(1,i);hits(1,i)+vside;hits(1,i)+vside;hits(1,i)];
        Ybound=[hits(2,i);hits(2,i);hits(2,i)+vside;hits(2,i)+vside;hits(2,i);hits(2,i);hits(2,i)+vside;hits(2,i)+vside];
        Zbound=[hits(3,i);hits(3,i);hits(3,i);hits(3,i);hits(3,i)+dz;hits(3,i)+dz;hits(3,i)+dz;hits(3,i)+dz];
        alpha=alphaShape(Xbound,Ybound,Zbound);
        
        color=Lij_slice(i)/max(Lij_slice)*[1,1,1];
        
        plot(alpha,'FaceColor',color)
        hold on
        
    end
end

%plot3([det(1,:),x],[det(2,:),y],[det(3,:),z],'r','LineWidth',5)
xlabel('Easting')
ylabel('Norhting')
zlabel('Elevation')