function makeTissue
%
%   Builds and saves a tissue in a rectangular cuboid voxel-space.
%   The tissue properties are loaded from makeTissueList.m.
%
%   First, define the wavelength in nm, and the tissue cuboid built of voxels.
%   Then, fill the voxel space with indices pointing to the required tissue
%   types as listed in makeTissueList.m.
%   This file produces the .mat file required as an input to
%   runMonteCarlo.m, and it displays both the tissue cuboid as well as an
%   overview over the optical, thermal and fluorescence properties of the tissue types (if available).
%   Note that the tissue cuboid is defined in xyz-coordinates (and not yxz).
%
%   Displays
%       Tissue cuboid
%       Tissue optical, thermal and fluorescence properties
%
%   Output
%       ./Data/[name].mat
%           file containing the 3D tissue cuboid definition (voxels)
%
%   Requires
%       deleteDataFiles.m
%       makeTissueList.m
%       lookMCmatlab.m
%

%% Acknowledgement
%   This function was inspired by maketissue of the mcxyz program hosted at omlc.org

%% Define parameters (user-specified)
wavelength  = 450;     % [nm] set the wavelength of the Monte Carlo simulation
name = 'tissue';        % name of the simulation
nx = 20;               % number of bins in the x direction
ny = 20;               % number of bins in the y direction
nz = 20;               % number of bins in the z direction
Lx = 1;                 % [cm] x size of simulation area
Ly = 1;                 % [cm] y size of simulation area
Lz = 1;                 % [cm] z size of simulation area

wavelength_fluorescence = NaN; % [nm] Fluorescence wavelength (set this to NaN for simulations without fluorescence)

%% Check for preexisting files
if(~deleteDataFiles(name)); return; end

%% Calculate x,y,z vectors and grids
dx = Lx/nx;             % [cm] size of x bins
dy = Ly/ny;             % [cm] size of y bins
dz = Lz/nz;             % [cm] size of z bins
x  = ((0:nx-1)-(nx-1)/2)*dx; % [cm] x position of centers of voxels
y  = ((0:ny-1)-(ny-1)/2)*dy; % [cm] y position of centers of voxels
z  = ((0:nz-1)+1/2)*dz;      % [cm] z position of centers of voxels
[X,Y,Z] = ndgrid(single(x),single(y),single(z)); % The single data type is used to conserve memory

%% Define tissue T(x,y,z) (user-specified)
%% Standard tissue test:
% T = 3*ones(nx,ny,nz,'uint8'); % "standard" tissue

%% Blood vessel example:
% zsurf = 0.01;
% epd_thick = 0.006;
% vesselradius  = 0.0100;
% vesseldepth = 0.04;
% T = 2*ones(nx,ny,nz,'uint8'); % fill background with water (gel)
% T(Z > zsurf) = 4; % epidermis
% T(Z > zsurf + epd_thick) = 5; % dermis
% T(X.^2 + (Z - (zsurf + vesseldepth)).^2 < vesselradius^2) = 6; % blood

%% Fluorescing cylinder example:
% cylinderradius  = 0.0100;
% T = 17*ones(nx,ny,nz,'uint8'); % fill background with fluorescence absorber
% T(Y.^2 + (Z - 3*cylinderradius).^2 < cylinderradius^2) = 16; % fluorescent tissue

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

%% Imaging example:
T = 1*ones(nx,ny,nz,'uint8'); % Air background
T(1:(nx*(ny+1)+1):end) = 18; % Set xyz diagonal positions to testscatterer
T(1:(nx*(ny+1)):end) = 18; % Set yz diagonal positions to testscatterer

%% Get the tissueList and the reduced T matrix
if(~isnan(wavelength_fluorescence))
    [~,tissueList_fluorescence] = makeTissueList(T,wavelength_fluorescence);
    [T, tissueList] = makeTissueList(T,wavelength);
    if(~any([tissueList.Y]>0))
        error('Fluorescence wavelength isn''t NaN, but none of the tissues have Y > 0');
    end
else
    [T, tissueList] = makeTissueList(T,wavelength);
end

%% Save output
if(~isnan(wavelength_fluorescence))
    save(['./Data/' name '.mat'],'dx','dy','dz','nx','ny','nz','x','y','z','tissueList','tissueList_fluorescence','T','wavelength','wavelength_fluorescence')
else
    save(['./Data/' name '.mat'],'dx','dy','dz','nx','ny','nz','x','y','z','tissueList','T','wavelength')
end
fprintf('./Data/%s.mat saved\n',name)

%% Make plots
lookMCmatlab(name);

end
