function defineGeometry(name)
%
%   Builds and saves a definition of the simulation geometry and the
%   optical media it contains in a rectangular cuboid voxel mesh.
%   The media properties are loaded from getMediaProperties.m.
%
%	Pay attention to the sections with headers that say "USER SPECIFIED:"
%	In those sections, you must fill in the parameters relevant for your simulation.
%	
%   Input
%       name
%           the basename of the file you want to store the geometry in
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
%       plotMCmatlab.m
%

%%%%%
%   Copyright 2017, 2018 by Dominik Marti and Anders K. Hansen, DTU Fotonik
%   This function was inspired by maketissue.m of the mcxyz program hosted at omlc.org
%
%   This file is part of MCmatlab.
%
%   MCmatlab is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
%
%   MCmatlab is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with MCmatlab.  If not, see <https://www.gnu.org/licenses/>.
%%%%%

%% USER SPECIFIED: Define simulation behavior
% Should the script run in silent mode? (disables overwrite prompt,
% command window text, progress indication and plot generation)
silentMode = false;

% Do you want to assume matched interfaces? If so, all refractive indices
% are assumed to be 1 and there is no Fresnel reflection or refraction.
% Otherwise, refractive indices from getMediaProperties are used. Note that
% non-matched interfaces must be normal to the z axis, so each xy-slice
% must have a constant refractive index. 
assumeMatchedInterfaces = true;

% Boundary type
% 0: No boundaries. Photons are allowed to leave the cuboid and are still
%    tracked outside, including absorption and scattering events. They get
%    terminated only if they wander too far (6 times the cuboid size).
% 1: Cuboid boundaries. All 6 cuboid surfaces are considered photon boundaries.
% 2: Top boundary only. Only the top surface (z = 0) is a photon boundary.
% Regardless of the boundary type, photons that wander 6 times the cuboid
% size will be terminated. When a photon hits a photon boundary at a position
% where the refractive index is 1, it escapes and may contribute to the
% signal of the light collector depending on its trajectory. Otherwise, the
% photon is just terminated, meaning that it cannot contribute to the light
% collector.
boundaryType = 1;


%% USER SPECIFIED: Define parameters
wavelength  = 532;		% [nm] set the wavelength of the Monte Carlo simulation
wavelength_f = NaN;		% [nm] Fluorescence wavelength (set this to NaN for simulations without fluorescence)

nx = 20;				% number of bins in the x direction
ny = 20;				% number of bins in the y direction
nz = 20;				% number of bins in the z direction
Lx = .1;				% [cm] x size of simulation area
Ly = .1;				% [cm] y size of simulation area
Lz = .1;				% [cm] z size of simulation area

%% Check for preexisting files
if(~silentMode && ~deleteDataFiles(name)); return; end

%% Calculate x,y,z vectors and grids
dx = Lx/nx;                  % [cm] size of x bins
dy = Ly/ny;                  % [cm] size of y bins
dz = Lz/nz;                  % [cm] size of z bins
x  = ((0:nx-1)-(nx-1)/2)*dx; % [cm] x position of centers of voxels
y  = ((0:ny-1)-(ny-1)/2)*dy; % [cm] y position of centers of voxels
z  = ((0:nz-1)+1/2)*dz;      % [cm] z position of centers of voxels
[X,Y,Z] = ndgrid(single(x),single(y),single(z)); % The single data type is used to conserve memory

%% USER SPECIFIED: Define medium matrix M(x,y,z)
% Fill the medium matrix M with indices that point to that location's
% medium type, as defined in getMediaProperties.
% Below are some examples you can get inspiration from.

%% Standard tissue example:
% M = 3*ones(nx,ny,nz,'uint8'); % "standard" tissue
% M(Z < 0.03) = 1; % Air

%% STL "dragon head" example
% M = ones(nx,ny,nz,'uint8');
% [STLvoxels,xSTL,zSTL,ySTL] = VOXELISE(nx,nz,ny,'Dragon_Head.stl'); % y and z dimensions intentionally swapped because we want to susequently rotate the object
% M(flip(flip(permute(STLvoxels,[1 3 2]),3),2)) = 3; % A few dimensions are rearranged to reorient the object and then "standard" tissue is assigned to the voxels inside the object
% 
% scalefactor = 1/500; % Factor to scale size of object by compared to the x, y, z data in the .stl file
% 
% % The voxel edge sizes and center positions have to be recalculated to fit the .stl voxelization
% dx = (xSTL(2)-xSTL(1))*scalefactor;
% dy = (ySTL(2)-ySTL(1))*scalefactor;
% dz = (zSTL(2)-zSTL(1))*scalefactor;
% x  = ((0:nx-1)-(nx-1)/2)*dx; % [cm] x position of centers of voxels
% y  = ((0:ny-1)-(ny-1)/2)*dy; % [cm] y position of centers of voxels
% z  = ((0:nz-1)+1/2)*dz;      % [cm] z position of centers of voxels

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
M = 1*ones(nx,ny,nz,'uint8'); % Air background
M(1:(nx*(ny+1)+1):end) = 18; % Set xyz diagonal positions to testscatterer
M(1:(nx*(ny+1)):end) = 18; % Set yz diagonal positions to testscatterer

%% Refraction and reflection example:
% M = ones(nx,ny,nz,'uint8'); % Air background
% M(Z>0.03) = 2; % Water
% M(Z>0.09) = 20; % Reflector

%% Get the mediaProperties and the reduced M matrix
if(~isnan(wavelength_f))
    [~,mediaProperties_f] = getMediaProperties(M,wavelength_f);
    [M, mediaProperties] = getMediaProperties(M,wavelength);
    if(~any([mediaProperties.Y]>0))
        error('Fluorescence wavelength isn''t NaN, but none of the media have Y > 0');
    end
else
    mediaProperties_f = NaN;
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
    [mediaProperties.n] = deal(1);
    RI = ones(nz,1);
end

%% Collect variables into a struct and save
G = struct('boundaryType',boundaryType,'dx',dx,'dy',dy,'dz',dz,'nx',nx,'ny',ny,'nz',nz,'x',x,'y',y,'z',z,...
    'wavelength',wavelength,'mediaProperties',mediaProperties,...
    'wavelength_f',wavelength_f,'mediaProperties_f',mediaProperties_f,...
    'M',M,'RI',RI);

save(['./Data/' name '.mat'],'G');
if(~silentMode); fprintf('./Data/%s.mat saved\n',name); end

%% Make plots
if(~silentMode); plotMCmatlab(name); end

end
