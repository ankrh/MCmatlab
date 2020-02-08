function model = defineGeometry(model)
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

%% Calculate basic cuboid variables and call geometry function to create the media matrix M
model.G.dx = model.G.Lx/model.G.nx;                  % [cm] size of x bins
model.G.dy = model.G.Ly/model.G.ny;                  % [cm] size of y bins
model.G.dz = model.G.Lz/model.G.nz;                  % [cm] size of z bins
model.G.x  = ((0:model.G.nx-1)-(model.G.nx-1)/2)*model.G.dx; % [cm] x position of centers of voxels
model.G.y  = ((0:model.G.ny-1)-(model.G.ny-1)/2)*model.G.dy; % [cm] y position of centers of voxels
model.G.z  = ((0:model.G.nz-1)+1/2)*model.G.dz;      % [cm] z position of centers of voxels
[X,Y,Z] = ndgrid(single(model.G.x),single(model.G.y),single(model.G.z)); % The single data type is used to conserve memory
model.G.M_raw = uint8(model.G.geomFunc(X,Y,Z,model.G.geomFuncParams));
end
