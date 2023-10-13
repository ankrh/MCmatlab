function model = plotMCmatlabHeat(model)
%   Requires
%       plotVolumetric.m
%
%   See also simulateHeatDistribution

%%%%%
%   Copyright 2017, 2018 by Dominik Marti and Anders K. Hansen, DTU Fotonik
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

G = model.G;
mP_fH = model.G.mediaPropertiesFunc(G.mediaPropParams); % Unsplit media
[unqMediaIdxs,~,M_trim] = unique(G.M_raw);
M_trim = reshape(M_trim,[G.nx G.ny G.nz]);
mP_fHtrim = mP_fH(unqMediaIdxs);
nM = numel(mP_fHtrim);
names = {mP_fHtrim.name};
numTemperatureSensors = size(model.HS.tempSensorPositions,1);

%% Write out the highest temperatures in each medium
fprintf('--------------------plotMCmatlabHeat---------------------\n');
for iM=1:nM
  if isfinite(model.HS.maxMediaTemps(iM))
    fprintf('Highest temperature obtained in %s is %.2f°C\n',names{iM},model.HS.maxMediaTemps(iM));
  end
end

%% Plot the geometry showing the temperature sensor locations and the sensor data
if numTemperatureSensors
  h_f = MCmatlab.NdimSliderPlot(G.M_raw,...
    'nFig',23,...
    'axisValues',{G.x,G.y,G.z},...
    'axisLabels',{'x [cm]','y [cm]','z [cm]','Geometry illustration'},...
    'indexLabels',names,...
    'linColormap',lines(256),...
    'axisEqual',true,...
    'reversedAxes',3,...
    'slicePositions',model.HS.slicePositions);
  h_f.Name = 'Temperature sensor illustration';
  title('Temperature sensor illustration');

  for i=numTemperatureSensors:-1:1
    indices = round((model.HS.tempSensorPositions(i,:)+[G.Lx G.Ly 0]/2)./[G.dx G.dy G.dz] + [0.5 0.5 0.5]); % +[1 1 1] to go from zero-index reference to one-index reference and -[0.5 0.5 0.5] to go from voxel corner position to voxel center positions
    indices = min([G.nx G.ny G.nz],max([1 1 1],indices)); % Coerce to the cuboid
    linindex = sub2ind(size(G.M_raw),indices(1),indices(2),indices(3));
    sensorNumbers{i,1} = num2str(i);
    sensorLabels{i,1} = [num2str(i) ', ' mP_fHtrim(G.M_raw(linindex)).name];
  end
  text(model.HS.tempSensorPositions(:,1),model.HS.tempSensorPositions(:,2),model.HS.tempSensorPositions(:,3),sensorNumbers,'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',18);

  h_f = figure(24);
  clf reset;
  set(h_f,'WindowStyle','Docked');
  h_f.Color = 'w';
  h_f.Name = 'Temperature sensor data';
  plot(model.HS.sensorsTimeVector,model.HS.sensorTemps,'LineWidth',2);
  set(gca,'FontSize',16);
  xlabel('Time [sec]');
  ylabel('Temperature [deg C]');
  title('Temperature sensor data');
  xlim(model.HS.sensorsTimeVector([1 end]));
  legend(sensorLabels,'Location','best');
  grid on;grid minor;
end

%% Plot thermal damage
% Let c be the concentration of undamaged molecules or living cells
% Omega = ln(c(t=0)/c(t=tau)) = int(A*exp(-E/(RT(t))),t=0..tau)
% c(0)/c(t) = exp(Omega) <=> c(t) = c(0)*exp(-Omega)
% So the fractional portion of damaged molecules/cells is:
% (c(0)-c(t))/c(0) = (c0-c0*exp(-Omega))/c0 = 1 - exp(-Omega)
if ~isnan(model.HS.Omega(1))
  M_damage = M_trim;
  M_damage(model.HS.Omega > 1) = nM + 1;
  namesWithDamage = [names 'damage'];
  h_f = MCmatlab.NdimSliderPlot(M_damage,...
    'nFig',25,...
    'axisValues',{G.x,G.y,G.z},...
    'axisLabels',{'x [cm]','y [cm]','z [cm]','Thermal damage illustration'},...
    'indexLabels',namesWithDamage,...
    'linColormap',lines(256),...
    'axisEqual',true,...
    'reversedAxes',3,...
    'slicePositions',model.HS.slicePositions);
  h_f.Name = 'Thermal damage illustration';
  title('Thermal damage illustration');
  fprintf('%.2e cm^3 was thermally damaged.\n',G.dx*G.dy*G.dz*sum(model.HS.Omega(:) > 1));
  
  h_f = MCmatlab.NdimSliderPlot(1-exp(-model.HS.Omega),...
    'nFig',26,...
    'axisValues',{G.x,G.y,G.z},...
    'axisLabels',{'x [cm]','y [cm]','z [cm]','Fractional damage'},...
    'axisEqual',true,...
    'reversedAxes',3,...
    'slicePositions',model.HS.slicePositions);
  h_f.Name = 'Fractional damage';
  title('Fractional damage (1 - e^{-\Omega})');
end
drawnow;
end

