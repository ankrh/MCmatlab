function plotMCmatlabHeat(HSinput,HSoutput)
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

if ~isfield(HSinput,'slicePositions')
  HSinput.slicePositions = [0.5 1 1];
end
if ~isfield(HSinput,'tempSensorPositions')
  HSinput.tempSensorPositions = [];
end

G = HSinput.G;
nM = length(G.mediaProperties); % Number of different media in simulation
numTemperatureSensors = size(HSinput.tempSensorPositions,1);

%% Plot the geometry showing the temperature sensor locations and the sensor data
if(numTemperatureSensors)
  geometryFigure = plotVolumetric(22,G.x,G.y,G.z,G.M,'MCmatlab_GeometryIllustration',G.mediaProperties,'slicePositions',HSinput.slicePositions);
  geometryFigure.Name = 'Temperature sensor illustration';
  title('Temperature sensor illustration');

  for i=numTemperatureSensors:-1:1
    indices = round((HSinput.tempSensorPositions(i,:)+[G.Lx G.Ly 0]/2)./[G.dx G.dy G.dz] + [0.5 0.5 0.5]); % +[1 1 1] to go from zero-index reference to one-index reference and -[0.5 0.5 0.5] to go from voxel corner position to voxel center positions
    indices = min([G.nx G.ny G.nz],max([1 1 1],indices)); % Coerce to the cuboid
    linindex = sub2ind(size(G.M),indices(1),indices(2),indices(3));
    sensorNumbers{i,1} = num2str(i);
    sensorLabels{i,1} = [num2str(i) ', ' G.mediaProperties(G.M(linindex)).name];
  end
  text(HSinput.tempSensorPositions(:,1),HSinput.tempSensorPositions(:,2),HSinput.tempSensorPositions(:,3),sensorNumbers,'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',18);

  if(~ishandle(23))
    temperatureSensorFigure = figure(23);
    temperatureSensorFigure.Position = [40 80 1100 650];
  else
    temperatureSensorFigure = figure(23);
  end
  temperatureSensorFigure.Color = 'w';
  clf;
  temperatureSensorFigure.Name = 'Temperature sensors';
  plot(HSoutput.sensorsTimeVector,HSoutput.sensorTemps,'LineWidth',2);
  set(gca,'FontSize',16);
  xlabel('Time [sec]');
  ylabel('Temperature [deg C]');
  title('Temperature sensors');
  xlim(HSoutput.sensorsTimeVector([1 end]));
  legend(sensorLabels,'Location','best');
  grid on;grid minor;
end

%% Plot thermal damage
if ~isnan(HSoutput.Omega(1))
  M_damage = G.M;
  M_damage(HSoutput.Omega > 1) = nM + 1;
  G.mediaProperties(nM + 1).name = 'damage';
  damageFigure = plotVolumetric(25,G.x,G.y,G.z,M_damage,'MCmatlab_GeometryIllustration',G.mediaProperties,'slicePositions',HSinput.slicePositions);
  damageFigure.Name = 'Thermal damage illustration';
  title('Thermal damage illustration');
  fprintf('%.2e cm^3 was thermally damaged.\n',G.dx*G.dy*G.dz*sum(sum(sum(HSoutput.Omega > 1))));
end
drawnow;
end

