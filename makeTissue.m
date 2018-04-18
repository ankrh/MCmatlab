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
%   Note that the tissue cube is defined in xyz-coordinates (and not yxz).
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
%   2018-04: Anders K. Hansen

%% Define parameters (user-specified)
wavelength  = 450;     % [nm] set the wavelength of the Monte Carlo simulation
name = 'tissue';        % name of the simulation
nx = 100;               % number of bins in the x direction
ny = 100;               % number of bins in the y direction
nz = 100;               % number of bins in the z direction
Lx = .1;                 % [cm] x size of simulation area
Ly = .1;                 % [cm] y size of simulation area
Lz = .1;                 % [cm] z size of simulation area

wavelength_fluorescence = 600; % [nm] Fluorescence wavelength (set this to NaN for simulations without fluorescence)

%% Calculate x,y,z vectors and grids
dx = Lx/nx;             % [cm] size of x bins
dy = Ly/ny;             % [cm] size of y bins
dz = Lz/nz;             % [cm] size of z bins
x  = ((0:nx-1)-(nx-1)/2)*dx;
y  = ((0:ny-1)-(ny-1)/2)*dy;
z  = ((0:nz-1)+1/2)*dz;
[X,Y,Z] = ndgrid(single(x),single(y),single(z)); % The single data type is used to conserve memory

%% Define tissue T(x,y,z) (user-specified)
%% Standard tissue test:
% T = 3*ones(nx,ny,nz,'uint8'); % "standard" tissue

%% Blood vessel example:
% zsurf = 0.01;
% epd_thick = 0.006;
% vesselradius  = 0.0100;
% T = ones(nx,ny,nz,'uint8'); % fill background with air
% T(Z > zsurf) = 4; % epidermis
% T(Z > zsurf + epd_thick) = 5; % dermis
% T(Y.^2 + (Z - (zsurf + 2*epd_thick + vesselradius)).^2 < vesselradius^2) = 6; % blood

%% Fluorescing cylinder example:
cylinderradius  = 0.0100;
T = 17*ones(nx,ny,nz,'uint8'); % fill background with fluorescent tissue
T(Y.^2 + (Z - 3*cylinderradius).^2 < cylinderradius^2) = 16; % fluorescence absorber

%% Hair example:
% zsurf = 0.02;  % position of gel/skin surface[cm]
% epd_thick = 0.01; % thickness of the epidermis [cm]
% hair_radius = 0.0075/2; % diameter varies from 17 - 180 micrometers, should increase with colouring and age
% hair_bulb_semiminor = 1.7*hair_radius; % [cm]
% hair_bulb_semimajor = sqrt(2)*hair_bulb_semiminor;
% hair_depth = 0.1; % varies from 0.06-0.3cm
% papilla_semiminor = hair_bulb_semiminor*5/12;
% papilla_semimajor = sqrt(2)*papilla_semiminor;
% 
% T = 2*ones(nx,ny,nz,'uint8'); % water (gel)
% T(Z > zsurf) = 4; % epidermis
% T(Z > zsurf+epd_thick) = 5; % dermis
% T(X.^2 + Y.^2 < hair_radius^2 & Z < zsurf+hair_depth) = 10; % hair
% T((X/hair_bulb_semiminor).^2 + (Y/hair_bulb_semiminor).^2 + ((Z-(zsurf+hair_depth))/hair_bulb_semimajor).^2 < 1) = 10; % hair
% T((X/papilla_semiminor).^2 + (Y/papilla_semiminor).^2 + ((Z-(zsurf+hair_depth+hair_bulb_semimajor-papilla_semimajor))/papilla_semimajor).^2 < 1) = 5; % dermis (papilla)

%% Solderpatch example:
% patch_radius        = 0.218;   	% [cm], cylinder radius
% patch_zi_start      = 1;
% patch_zi_end        = 5;
% vessel_radius       = 0.19;   	% [cm], cylinder radius
% water_radius        = 0.15;   	% [cm], cylinder radius
% fibre_radius        = 0.04;   	% [cm], cylinder radius
% 
% T = ones(nx,ny,nz,'uint8'); % fill background with air
% T(X.^2 + Y.^2 < patch_radius^2 & Z >= patch_zi_start & Z <= patch_zi_end) = 12; % patch
% T(X.^2 + Y.^2 < vessel_radius^2) = 7; % vessel
% T(X.^2 + Y.^2 < water_radius^2) = 2; % water
% T(X.^2 + Y.^2 < fibre_radius^2) = 11; % fibre

%% Get the tissueList and reduce T
if(~isnan(wavelength_fluorescence))
    [~,tissueList_fluorescence] = makeTissueList(T,wavelength_fluorescence);
    [T, tissueList] = makeTissueList(T,wavelength);
    if(~any([tissueList.Y]>0))
        error('Fluorescence wavelength isn''t NaN, but none of the tissues have Y > 0');
    end
else
    [T, tissueList] = makeTissueList(T,wavelength);
end

%% Check for preexisting files
if(exist(['./Data/' name '.mat'],'file') || exist(['./Data/' name '_MCoutput.mat'],'file') || exist(['./Data/' name '_MCoutput_fluorescence.mat'],'file') || exist(['./Data/' name '_heatSimoutput.mat'],'file') || exist(['./Data/' name '_heatSimoutput.mp4'],'file'))
    if(strcmp(questdlg('Tissue definition and/or computation results by this name already exist. Delete existing files?','Overwrite prompt','Yes','No, abort','Yes'),'No, abort'))
        fprintf('Aborted without saving data.\n');
        return;
    end
    
    if(exist(['./Data/' name '_MCoutput.mat'],'file'))
        delete(['./Data/' name '_MCoutput.mat']);
    end
    if(exist(['./Data/' name '_MCoutput_fluorescence.mat'],'file'))
        delete(['./Data/' name '_MCoutput_fluorescence.mat']);
    end
    if(exist(['./Data/' name '_heatSimoutput.mat'],'file'))
        delete(['./Data/' name '_heatSimoutput.mat']);
    end
    if(exist(['./Data/' name '_heatSimoutput.mp4'],'file'))
        delete(['./Data/' name '_heatSimoutput.mp4']);
    end
end

%% Save output
if(~isnan(wavelength_fluorescence))
    save(['./Data/' name '.mat'],'x','y','z','tissueList','tissueList_fluorescence','T','wavelength','wavelength_fluorescence')
else
    save(['./Data/' name '.mat'],'x','y','z','tissueList','T','wavelength')
end
fprintf('./Data/%s.mat saved\n',name)

%% Make plots
lookmcxyz(name);

return
