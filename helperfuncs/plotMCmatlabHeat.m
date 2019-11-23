function plotMCmatlabHeat(model)
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

G = model.G;
mP_fH = model.HS.mediaProperties_funcHandles;
nM = length(mP_fH); % Number of different media in simulation
numTemperatureSensors = size(model.HS.tempSensorPositions,1);

%% Show thermal media properties
h_f = plotMediaProperties(22,model,3);
h_f.Name = 'Thermal media properties';

%% Write out the highest temperatures in each medium
fprintf('--------------------plotMCmatlabHeat---------------------\n');
for idx=1:nM
  if isfinite(model.HS.maxMediaTemps(idx))
    fprintf('Highest temperature obtained in %s is %.2f°C\n',mP_fH(idx).name,model.HS.maxMediaTemps(idx));
  end
end

%% Plot the geometry showing the temperature sensor locations and the sensor data
if numTemperatureSensors
  geometryFigure = plotVolumetric(23,G.x,G.y,G.z,G.M_raw,'MCmatlab_GeometryIllustration',mP_fH,'slicePositions',model.HS.slicePositions);
  geometryFigure.Name = 'Temperature sensor illustration';
  title('Temperature sensor illustration');

  for i=numTemperatureSensors:-1:1
    indices = round((model.HS.tempSensorPositions(i,:)+[G.Lx G.Ly 0]/2)./[G.dx G.dy G.dz] + [0.5 0.5 0.5]); % +[1 1 1] to go from zero-index reference to one-index reference and -[0.5 0.5 0.5] to go from voxel corner position to voxel center positions
    indices = min([G.nx G.ny G.nz],max([1 1 1],indices)); % Coerce to the cuboid
    linindex = sub2ind(size(G.M_raw),indices(1),indices(2),indices(3));
    sensorNumbers{i,1} = num2str(i);
    sensorLabels{i,1} = [num2str(i) ', ' mP_fH(G.M_raw(linindex)).name];
  end
  text(model.HS.tempSensorPositions(:,1),model.HS.tempSensorPositions(:,2),model.HS.tempSensorPositions(:,3),sensorNumbers,'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',18);

  if ~ishandle(24)
    temperatureSensorFigure = figure(24);
    temperatureSensorFigure.Position = [40 160 1100 650];
  else
    temperatureSensorFigure = figure(24);
  end
  temperatureSensorFigure.Color = 'w';
  clf;
  temperatureSensorFigure.Name = 'Temperature sensor data';
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
  M_damage = G.M_raw;
  M_damage(model.HS.Omega > 1) = nM + 1;
  mP_fH(nM + 1).name = 'damage';
  damageFigure = plotVolumetric(25,G.x,G.y,G.z,M_damage,'MCmatlab_GeometryIllustration',mP_fH,'slicePositions',model.HS.slicePositions);
  damageFigure.Name = 'Thermal damage illustration';
  title('Thermal damage illustration');
  fprintf('%.2e cm^3 was thermally damaged.\n',G.dx*G.dy*G.dz*sum(sum(sum(model.HS.Omega > 1))));
  
  FDfigure = plotVolumetric(26,G.x,G.y,G.z,1-exp(-model.HS.Omega),'MCmatlab_fromZero','slicePositions',model.HS.slicePositions);
  FDfigure.Name = 'Fractional damage';
  title('Fractional damage (1 - e^{-\Omega})');
end
drawnow;
end

