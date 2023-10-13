function model = plotMCmatlabGeom(model)
%   Displays
%       Geometry cuboid
%       Media optical, thermal and fluorescence properties
%
%   Requires
%       NdimSliderPlot.m
%

%%%%%
%   Copyright 2017, 2018 by Dominik Marti and Anders K. Hansen, DTU Fotonik
%   This function was inspired by lookmcxyz.m of the mcxyz MC program hosted at omlc.org
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

if strcmp(matlabdrive,'/MATLAB Drive') && strcmp(matlabroot, '/MATLAB')
  % dirty hack to check whether we are on MATLAB online, 
  % where the Figures window tabs can't be on the side of the figure window
  com.mathworks.mde.desk.MLDesktop.getInstance.setDocumentBarPosition('Figures',1); % Set Figures window tabs to be on top
else
  com.mathworks.mde.desk.MLDesktop.getInstance.setDocumentBarPosition('Figures',7); % Set Figures window tabs to be on left side
end
model.G = model.G.updateGeometry;

%% Make geometry plot
mediaProperties = model.G.mediaPropertiesFunc(model.G.mediaPropParams);
if ~isa(mediaProperties,'MCmatlab.mediumProperties')
  error('Error: The syntax for defining media properties has changed slightly. You must now put "mediaProperties = MCmatlab.mediumProperties;" as the first line of the mediaProperties function. See the examples for details.');
end
[uniqueMedia,~,M_trimmed] = unique(model.G.M_raw);
M_trimmed = reshape(M_trimmed,model.G.nx,model.G.ny,model.G.nz);
mediaProperties_trimmed = mediaProperties(uniqueMedia);
h_f = MCmatlab.NdimSliderPlot(M_trimmed,...
  'nFig',1,...
  'axisValues',{model.G.x,model.G.y,model.G.z},...
  'axisLabels',{'x [cm]','y [cm]','z [cm]','Geometry illustration'},...
  'indexLabels',{mediaProperties_trimmed.name},...
  'linColormap',lines(256),...
  'axisEqual',true,...
  'reversedAxes',3);
h_f.Name = 'Geometry illustration';

drawnow;
figure(1);
end
