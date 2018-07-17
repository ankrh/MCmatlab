function G = defineGeometry(name,varargin)
%
%   Builds and saves a definition of the simulation geometry and the
%   optical media it contains in a rectangular cuboid voxel-space.
%   The media properties are loaded from getMediaProperties.m.
%
%   First, define the wavelength in nm, and the geometry cuboid built of voxels.
%   Then, fill the voxel space with indices pointing to the required media
%   as listed in getMediaProperties.m.
%   This file produces the .mat file required as an input to
%   runMonteCarlo.m, and it displays both the geometry cuboid as well as an
%   overview over the optical, thermal and fluorescence properties of the media (if available).
%   Note that the geometry cuboid is defined in xyz-coordinates (and not yxz).
%
%   Input
%       name
%           the basename of the file you want to store the geometry in
%       varargin
%           if 'silent' is specified as an additional argument, disables
%           overwrite prompt and command window text
%
%   Displays
%       Geometry cuboid
%       Media optical, thermal and fluorescence properties
%
%   Output
%       ./Data/[name].mat
%           file containing the 3D geometry cuboid definition (voxels)
%
%   Requires
%       deleteDataFiles.m
%       getMediaProperties.m
%       lookMCmatlab.m
%

%% Acknowledgement
%   This function was inspired by maketissue of the mcxyz program hosted at omlc.org

%% USER SPECIFIED: Define parameters
wavelength  = 450;     % [nm] set the wavelength of the Monte Carlo simulation
wavelength_fluorescence = NaN; % [nm] Fluorescence wavelength (set this to NaN for simulations without fluorescence)

nx = 100;               % number of bins in the x direction
ny = 100;               % number of bins in the y direction
nz = 100;               % number of bins in the z direction
Lx = .1;                 % [cm] x size of simulation area
Ly = .1;                 % [cm] y size of simulation area
Lz = .1;                 % [cm] z size of simulation area

% Do you want to assume matched interfaces? If so, there is no Fresnel
% reflection or refraction. Otherwise, refractive indices from
% getMediaProperties are used. Note that non-matched interfaces must be normal
% to the z axis, so each xy-slice must have a constant refractive index.
assumeMatchedInterfaces = false; 

% Boundary type
% 0: No boundaries. Photons wander freely also outside the cuboid and
%    get killed only if they wander too far (6 times the cuboid size).
% 1: Escape at boundaries. Photons that stray outside the cuboid get
%    killed immediately.
% 2: Escape at surface only. Photons that hit the top surface get killed
%    immediately, photons hitting other surfaces can wander up to 6 times
%    the cuboid size.
boundaryFlag = 1;

%% Check if silent mode was specified
silentMode = any(strcmpi(varargin,'silent'));

%% Check for preexisting files
if(~silentMode && ~deleteDataFiles(name)); return; end

%% Calculate x,y,z vectors and grids
dx = Lx/nx;             % [cm] size of x bins
dy = Ly/ny;             % [cm] size of y bins
dz = Lz/nz;             % [cm] size of z bins
x  = ((0:nx-1)-(nx-1)/2)*dx; % [cm] x position of centers of voxels
y  = ((0:ny-1)-(ny-1)/2)*dy; % [cm] y position of centers of voxels
z  = ((0:nz-1)+1/2)*dz;      % [cm] z position of centers of voxels
[X,Y,Z] = ndgrid(single(x),single(y),single(z)); % The single data type is used to conserve memory

%% USER SPECIFIED: Define medium matrix M(x,y,z)
%% Standard tissue test:
% M = 3*ones(nx,ny,nz,'uint8'); % "standard" tissue

%% Blood vessel example:
% zsurf = 0.01;
% epd_thick = 0.006;
% vesselradius  = 0.0100;
% vesseldepth = 0.04;
% M = 2*ones(nx,ny,nz,'uint8'); % fill background with water (gel)
% M(Z > zsurf) = 4; % epidermis
% M(Z > zsurf + epd_thick) = 5; % dermis
% M(X.^2 + (Z - (zsurf + vesseldepth)).^2 < vesselradius^2) = 6; % blood

%% Fluorescing cylinder example:
% cylinderradius  = 0.0100;
% M = 17*ones(nx,ny,nz,'uint8'); % fill background with fluorescence absorber
% M(Y.^2 + (Z - 3*cylinderradius).^2 < cylinderradius^2) = 16; % fluorescer

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
% M = 2*ones(nx,ny,nz,'uint8'); % water (gel)
% M(Z > zsurf) = 4; % epidermis
% M(Z > zsurf+epd_thick) = 5; % dermis
% M(X.^2 + Y.^2 < hair_radius^2 & Z < zsurf+hair_depth) = 10; % hair
% M((X/hair_bulb_semiminor).^2 + (Y/hair_bulb_semiminor).^2 + ((Z-(zsurf+hair_depth))/hair_bulb_semimajor).^2 < 1) = 10; % hair
% M((X/papilla_semiminor).^2 + (Y/papilla_semiminor).^2 + ((Z-(zsurf+hair_depth+hair_bulb_semimajor-papilla_semimajor))/papilla_semimajor).^2 < 1) = 5; % dermis (papilla)

%% Solderpatch example:
% patch_radius        = 0.218;   	% [cm], cylinder radius
% patch_zi_start      = 1;
% patch_zi_end        = 5;
% vessel_radius       = 0.19;   	% [cm], cylinder radius
% water_radius        = 0.15;   	% [cm], cylinder radius
% fibre_radius        = 0.04;   	% [cm], cylinder radius
% 
% M = ones(nx,ny,nz,'uint8'); % fill background with air
% M(X.^2 + Y.^2 < patch_radius^2 & Z >= patch_zi_start & Z <= patch_zi_end) = 12; % patch
% M(X.^2 + Y.^2 < vessel_radius^2) = 7; % vessel
% M(X.^2 + Y.^2 < water_radius^2) = 2; % water
% M(X.^2 + Y.^2 < fibre_radius^2) = 11; % fibre

%% Imaging example:
% M = 1*ones(nx,ny,nz,'uint8'); % Air background
% M(1:(nx*(ny+1)+1):end) = 18; % Set xyz diagonal positions to testscatterer
% M(1:(nx*(ny+1)):end) = 18; % Set yz diagonal positions to testscatterer

%% Refraction and reflection example:
M = ones(nx,ny,nz,'uint8'); % Air background
M(Z>0.03) = 2; % Water
M(Z>0.09) = 20; % Reflector

%% Get the mediaProperties and the reduced M matrix
if(~isnan(wavelength_fluorescence))
    [~,mediaProperties_fluorescence] = getMediaProperties(M,wavelength_fluorescence);
    [M, mediaProperties] = getMediaProperties(M,wavelength);
    if(~any([mediaProperties.Y]>0))
        error('Fluorescence wavelength isn''t NaN, but none of the media have Y > 0');
    end
else
    mediaProperties_fluorescence = NaN;
    [M, mediaProperties] = getMediaProperties(M,wavelength);
end

%% Extract the refractive indices 
if(~assumeMatchedInterfaces)
    for j=1:length(mediaProperties) % Check that all media have a refractive index defined
        if(~isfield(mediaProperties,'n') || any(isempty(mediaProperties(j).n)))
            error('assumeMatchedInterfaces is false, but refractive index isn''t defined for all media');
        end
    end
    n_vec = [mediaProperties.n];
    for j=1:nz % Check that each xy slice has constant refractive index, so refractive index is only a function of z
        if(length(unique(n_vec(M(:,:,j)))) > 1)
            error('assumeMatchedInterfaces is false, but refractive index isn''t constant for z index %d (z = %f).\nEach xy slice must have constant refractive index.',j,z(j));
        end
    end
    RI = n_vec(M(1,1,:));
else
    RI = NaN;
end

%% Collect variables into a struct and save
G = struct('boundaryFlag',boundaryFlag,'dx',dx,'dy',dy,'dz',dz,'nx',nx,'ny',ny,'nz',nz,'x',x,'y',y,'z',z,...
    'wavelength',wavelength,'mediaProperties',mediaProperties,...
    'wavelength_fluorescence',wavelength_fluorescence,'mediaProperties_fluorescence',mediaProperties_fluorescence,...
    'M',M,'RI',RI);

save(['./Data/' name '.mat'],'G');
if(~silentMode) fprintf('./Data/%s.mat saved\n',name); end

%% Make plots
if(~silentMode) lookMCmatlab(name); end

end
