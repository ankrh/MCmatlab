function makeTissue
%
%   Builds and saves a tissue in a voxel-space. The tissue properties are
%   loaded from makeTissueList.m.
%
%   First, define the wavelength in nm, and the tissue cube built of voxels.
%   Then, fill the voxel space with indeces pointing to the required tissue
%   types as listed in makeTissueList.m.
%   This file produces the .mat file required as an input to
%   runMonteCarlo.m, and it displays both the tissue cube as well as an
%   overview over the optical and thermal properties of the tissue types.
%
%   Displays
%       Tissue cube
%       Tissue optical and thermal properties
%
%   Output
%       ./Data/[name].mat
%           file containing the 3D tissue cube definition (voxels)
%
%   Requires
%       makeTissueList.m
%       plotVolumetric.m
%       plotTissueProperties.m
%

%% Updates
%   2014-08: Steven L. Jacques
%   2017-06: Anders K. Hansen & Dominik Marti, DTU Fotonik

%% Define parameters (user-specified)
wavelength  = 1000;     % [nm] set the wavelength of the Monte Carlo simulation
name = 'dentin';        % name of the simulation
nx = 100;               % number of bins in the x direction
ny = 100;               % number of bins in the y direction
nz = 100;               % number of bins in the z direction
Lx = 1;                 % [cm] x size of simulation area
Ly = 1;                 % [cm] y size of simulation area
Lz = 1;                 % [cm] z size of simulation area

%% Calculate x,y,z vectors and grids
dx = Lx/nx;             % [cm] size of x bins
dy = Ly/ny;             % [cm] size of y bins
dz = Lz/nz;             % [cm] size of z bins
x  = ((0:nx-1)-(nx-1)/2)*dx;
y  = ((0:ny-1)-(ny-1)/2)*dy;
z  = ((0:nz-1)+1/2)*dz;
[X,Y,Z] = ndgrid(single(x),single(y),single(z)); % The single data type is used to conserve memory

%% Define tissue T(x,y,z) (user-specified)
% Dentin example:
T = 3*ones(nx,ny,nz,'uint8'); % fill background with air
T(Y.^2 + (Z-nz*dz/2).^2 < 0.300^2) = 2; % enamel
T(Y.^2 + (Z-nz*dz/2).^2 < 0.170^2) = 1; % dentin
T(Y.^2 + (Z-nz*dz/2).^2 < 0.050^2) = 5; % blood

% % Solderpatch example:
% patch_radius        = 0.218;   	% [cm], cylinder radius
% patch_zi_start      = 1;
% patch_zi_end        = 5;
% vessel_radius       = 0.19;   	% [cm], cylinder radius
% water_radius        = 0.15;   	% [cm], cylinder radius
% fibre_radius        = 0.04;   	% [cm], cylinder radius
% 
% T = 3*ones(nx,ny,nz,'uint8'); % fill background with air
% T(X.^2 + Y.^2 < patch_radius^2 & Z >= patch_zi_start & Z <= patch_zi_end) = 8; % patch
% T(X.^2 + Y.^2 < vessel_radius^2) = 7; % vessel
% T(X.^2 + Y.^2 < water_radius^2) = 4; % water
% T(X.^2 + Y.^2 < fibre_radius^2) = 6; % fibre

% Hair example:
% zsurf = 0.02;  % position of gel/skin surface[cm]
% epd_thick = 0.01; % thickness of the epidermis [cm]
% hair_radius = 0.0075/2; % diameter varies from 17 - 180 micrometers, should increase with colouring and age
% hair_bulb_semiminor = 1.7*hair_radius; % [cm]
% hair_bulb_semimajor = sqrt(2)*hair_bulb_semiminor;
% hair_depth = 0.1; % varies from 0.06-0.3cm
% papilla_semiminor = hair_bulb_semiminor*5/12;
% papilla_semimajor = sqrt(2)*papilla_semiminor;
% 
% T = 4*ones(nx,ny,nz,'uint8'); % water (gel)
% T(Z > zsurf) = 10; % epidermis
% T(Z > zsurf+epd_thick) = 9; % dermis
% T(X.^2 + Y.^2 < hair_radius^2 & Z < zsurf+hair_depth) = 15; % hair
% T((X/hair_bulb_semiminor).^2 + (Y/hair_bulb_semiminor).^2 + ((Z-(zsurf+hair_depth))/hair_bulb_semimajor).^2 < 1) = 15; % hair
% T((X/papilla_semiminor).^2 + (Y/papilla_semiminor).^2 + ((Z-(zsurf+hair_depth+hair_bulb_semimajor-papilla_semimajor))/papilla_semimajor).^2 < 1) = 9; % dermis (papilla)

% Blood vessel example:
% zsurf = 0.01;
% epd_thick = 0.006;
% vesselradius  = 0.0100;
% T = 3*ones(nx,ny,nz,'uint8'); % fill background with air
% T(Z > zsurf) = 10; % epidermis
% T(Z > zsurf + epd_thick) = 9; % dermis
% T(Y.^2 + (Z - (zsurf + 2*epd_thick + vesselradius)).^2 < vesselradius^2) = 5; % blood

% Standard tissue test:
% T = 14*ones(nx,ny,nz,'uint8'); % "standard" tissue

%% Discard the unused tissue types
tissueList = makeTissueList(wavelength);
nT = length(unique(T)); % Number of different tissues in simulation
tissueMap = zeros(1,length(tissueList),'uint8');
tissueMap(unique(T)) = 1:nT;
tissueList = tissueList(unique(T));
T = tissueMap(T); % Reduced tissue matrix, using only numbers from 1 up to the number of tissues actually used

%% Make plots
figure(1);clf;
plotVolumetric(x,y,z,T,tissueList)
title('Tissue type illustration');
drawnow;

figure(7);clf;
plotTissueProperties(tissueList);
drawnow;

%% Save output
save(['./Data/' name '.mat'],'x','y','z','tissueList','T','wavelength')
fprintf('./Data/%s.mat saved\n',name)

return
