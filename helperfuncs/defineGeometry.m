function G = defineGeometry(G)
%   Created 2018 by Dominik Marti and Anders K. Hansen, DTU Fotonik
%   
%   This function was inspired by maketissue.m of the mcxyz program hosted at omlc.org
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

if ~isfield(G,'wavelength_f')
	G.wavelength_f = NaN;
end
if ~isfield(G,'GeomFuncParams')
	G.GeomFuncParams = {};
end

G.dx = G.Lx/G.nx;                  % [cm] size of x bins
G.dy = G.Ly/G.ny;                  % [cm] size of y bins
G.dz = G.Lz/G.nz;                  % [cm] size of z bins
G.x  = ((0:G.nx-1)-(G.nx-1)/2)*G.dx; % [cm] x position of centers of voxels
G.y  = ((0:G.ny-1)-(G.ny-1)/2)*G.dy; % [cm] y position of centers of voxels
G.z  = ((0:G.nz-1)+1/2)*G.dz;      % [cm] z position of centers of voxels
[X,Y,Z] = ndgrid(single(G.x),single(G.y),single(G.z)); % The single data type is used to conserve memory
G.M = uint8(G.GeomFunc(X,Y,Z,G.GeomFuncParams));

%% Get the mediaProperties and the reduced M matrix
if(~isnan(G.wavelength_f))
    [~,G.mediaProperties_f] = getMediaProperties(G.M,G.wavelength_f);
    [G.M, G.mediaProperties] = getMediaProperties(G.M,G.wavelength);
    if(~any([G.mediaProperties.Y]>0))
        error('Fluorescence wavelength isn''t NaN, but none of the media have Y > 0');
    end
else
    G.mediaProperties_f = NaN;
    [G.M, G.mediaProperties] = getMediaProperties(G.M,G.wavelength);
end

%% Extract the refractive indices 
if(~G.matchedInterfaces)
    for j=1:length(G.mediaProperties) % Check that all media have a refractive index defined
        if(~isfield(G.mediaProperties,'n') || any(isempty(G.mediaProperties(j).n)))
            error('matchedInterfaces is false, but refractive index isn''t defined for all media');
        end
    end
    n_vec = [G.mediaProperties.n];
    for j=1:G.nz % Check that each xy slice has constant refractive index, so refractive index is only a function of z
        if(length(unique(n_vec(G.M(:,:,j)))) > 1)
            error('matchedInterfaces is false, but refractive index isn''t constant for z index %d (z = %f).\nEach xy slice must have constant refractive index.',j,z(j));
        end
    end
    G.RI = n_vec(G.M(1,1,:));
else
    [G.mediaProperties.n] = deal(1);
    G.RI = ones(G.nz,1);
end

if(~G.silentMode); plotMCmatlabGeom(G); end
end
