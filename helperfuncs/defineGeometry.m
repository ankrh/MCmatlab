function G = defineGeometry(G)
%   Builds and saves a definition of the simulation geometry and the
%   optical media it contains in a rectangular cuboid voxel mesh.
%   The media properties are loaded from getMediaProperties.m.
%
%   Requires
%       getMediaProperties.m
%
%   See also getMediaProperties, mediaPropertiesLibrary, runMonteCarlo

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

%% Use default values for unspecified fields
if ~isfield(G,'silentMode')
  G.silentMode = false;
end
if ~isfield(G,'GeomFuncParams')
  G.GeomFuncParams = {};
end
if ~isfield(G,'mediaPropParams')
  G.mediaPropParams = {};
end

%% Calculate basic cuboid variables and call geometry function to create the media matrix M
G.dx = G.Lx/G.nx;                  % [cm] size of x bins
G.dy = G.Ly/G.ny;                  % [cm] size of y bins
G.dz = G.Lz/G.nz;                  % [cm] size of z bins
G.x  = ((0:G.nx-1)-(G.nx-1)/2)*G.dx; % [cm] x position of centers of voxels
G.y  = ((0:G.ny-1)-(G.ny-1)/2)*G.dy; % [cm] y position of centers of voxels
G.z  = ((0:G.nz-1)+1/2)*G.dz;      % [cm] z position of centers of voxels
[X,Y,Z] = ndgrid(single(G.x),single(G.y),single(G.z)); % The single data type is used to conserve memory
G.M = uint8(G.GeomFunc(X,Y,Z,G.GeomFuncParams));
end
