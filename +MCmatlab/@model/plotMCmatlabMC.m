function model = plotMCmatlabMC(model,varargin)
%   Displays (if calculated)
%       Absorbed power
%       Fluence rate of all photons
%       Example paths of some of the photons
%       Distribution of escaped photons in the far field
%       Fluence rates (intensities) on all boundaries
%   And, if a light collector was defined, displays (if calculated)
%       An illustration of the light collector angle and placement
%       Image generated (which might be time-resolved)
%       Fluence rate of photons that hit the light collector
%
%   Requires
%       NdimSliderPlot.m
%       NdimSliderPlotRedraw.m
%
%   See also runMonteCarlo

%%%%%
%   Copyright 2017, 2018 by Dominik Marti and Anders K. Hansen
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

% 1: Geometry illustration
% 2: Media properties
% 3: Absorption
% 4: NFR
% 5: Paths
% 6: LC illustration
% 7: NFRdet
% 8: Collected light
% 9: Far field
% 10: Boundary Irradiances
% 11: Emitter distribution
% 12-20: Same as incident light but for fluorescence
% 21: Temperature plot during heat sim
% 22: Thermal media properties
% 23: Temperature sensor locations
% 24: Temperature sensor data
% 25: Thermal damage illustration
% 26: Fractional damage plot
% 31: STL shape

com.mathworks.mde.desk.MLDesktop.getInstance.setDocumentBarPosition('Figures',7); % Set Figures window tabs to be on left side
model.G = model.G.updateGeometry;

G = model.G;

simFluorescence = ~isempty(varargin) && strcmp(varargin{1},'fluorescence');

if simFluorescence
  simType = 2;
  MCorFMC = model.FMC;
  figNumOffset = 10;
  fluorescenceOrIncident = 'fluorescence ';
  fluorescenceOrNothing = 'fluorescence ';
  fprintf('----------------Fluorescence plotMCmatlab-----------------\n');
else
  simType = 1;
  MCorFMC = model.MC;
  figNumOffset = 0;
  fluorescenceOrIncident = 'incident ';
  fluorescenceOrNothing = '';
  fprintf('----------------------plotMCmatlab-----------------------\n');
end

if isnan(MCorFMC.nPhotons)
  error('Error: The simulation has not yet been run.');
end
if MCorFMC.nPhotons == 0
  error('Error: No photons were successfully simulated. Check your model file to ensure photons launch within the simulation cuboid.');
end

LC = MCorFMC.LC;

%% Plot emitter distribution
P_in = 1;
if simFluorescence
  h_f = MCmatlab.NdimSliderPlot(model.FMC.sourceDistribution,...
    'nFig',11,...
    'axisValues',{G.x,G.y,G.z,MCorFMC.wavelength},...
    'axisLabels',{'x [cm]','y [cm]','z [cm]','λ [nm]','Fluorescence emitter distribution [W/cm^3/W.incident]'},...
    'fromZero',true,...
    'axisEqual',true,...
    'reversedAxes',3);
  set(h_f,'WindowStyle','Docked');
  h_f.Name = 'Fluorescence emitters';

  P_exc_abs = G.dx*G.dy*G.dz*sum(model.MC.NA(:));
  P_flu_emit = G.dx*G.dy*G.dz*sum(model.FMC.sourceDistribution(:));
  fprintf('%.3g%% of absorbed incident light was re-emitted as fluorescence.\n',100*P_flu_emit/P_exc_abs);
  P_in = P_flu_emit;
end

%% Make media properties plot
h_f = plotMediaProperties(2,model);
set(h_f,'WindowStyle','Docked');
h_f.Name = 'Media properties';

if ~isscalar(MCorFMC.NFR)
  %% Make absorption plot
  h_f = MCmatlab.NdimSliderPlot(MCorFMC.NA,...
    'nFig',3 + figNumOffset,...
    'axisValues',{G.x,G.y,G.z,MCorFMC.wavelength},...
    'axisLabels',{'x [cm]','y [cm]','z [cm]','λ [nm]',['Normalized ' fluorescenceOrNothing 'absorbed power per unit volume [W/cm^3/W.incident]']},...
    'fromZero',true,...
    'axisEqual',true,...
    'reversedAxes',3);
  set(h_f,'WindowStyle','Docked');
  h_f.Name = ['Normalized ' fluorescenceOrNothing 'absorption'];
  
  fprintf(['%.3g%% of ' fluorescenceOrIncident 'light was absorbed within the cuboid.\n'],100*G.dx*G.dy*G.dz*sum(MCorFMC.NA(:))/P_in);
  
  %% Make fluence rate plot
  h_f = MCmatlab.NdimSliderPlot(MCorFMC.NFR,...
    'nFig',4 + figNumOffset,...
    'axisValues',{G.x,G.y,G.z,MCorFMC.wavelength},...
    'axisLabels',{'x [cm]','y [cm]','z [cm]','λ [nm]',['Normalized ' fluorescenceOrNothing 'fluence rate [W/cm^2/W.incident]']},...
    'fromZero',true,...
    'axisEqual',true,...
    'reversedAxes',3);
  set(h_f,'WindowStyle','Docked');
  h_f.Name = ['Normalized ' fluorescenceOrNothing 'fluence rate'];
end

%% Plot example paths
if MCorFMC.nExamplePaths > 0
  h_f = MCmatlab.NdimSliderPlot(G.M_raw,...
    'nFig',5 + figNumOffset,...
    'axisValues',{G.x,G.y,G.z},...
    'axisLabels',{'x [cm]','y [cm]','z [cm]','Geometry Illustration'},...
    'indexLabels',MCorFMC.mediaProperties.name,...
    'linColormap',lines(256),...
    'axisEqual',true,...
    'reversedAxes',3);
  legend off;
  set(h_f,'WindowStyle','Docked');
  if simFluorescence
    h_f.Name = 'Fluorescence photon paths';
  else
    h_f.Name = 'Photon paths';
  end
  title(h_f.Name);
  box on;grid on;grid minor;

  previousNaNidx = 1;
  linenumber = 0;
  for idx=3:(size(MCorFMC.examplePaths,2)+1)
    if idx == size(MCorFMC.examplePaths,2)+1 || isnan(MCorFMC.examplePaths(1,idx))
      linenumber = linenumber + 1;
      xdata = MCorFMC.examplePaths(1,previousNaNidx+1:idx-1);
      ydata = MCorFMC.examplePaths(2,previousNaNidx+1:idx-1);
      zdata = MCorFMC.examplePaths(3,previousNaNidx+1:idx-1);
      adata = MCorFMC.examplePaths(4,previousNaNidx+1:idx-1);
      previousNaNidx = idx;
      surface([xdata;xdata],...
              [ydata;ydata],...
              [zdata;zdata],...
              'AlphaData',uint8(64*[adata;adata]),...
              'EdgeColor',[0 0 0],...
              'EdgeAlpha','flat',...
              'FaceColor','none',...
              'LineWidth',2);
    end
  end
end

%% If there's a light collector, show its orientation and the detected light
if MCorFMC.useLightCollector
  h_f = MCmatlab.NdimSliderPlot(G.M_raw,...
    'nFig',6 + figNumOffset,...
    'axisValues',{G.x,G.y,G.z},...
    'axisLabels',{'x [cm]','y [cm]','z [cm]','Geometry Illustration'},...
    'indexLabels',MCorFMC.mediaProperties.name,...
    'linColormap',lines(256),...
    'axisEqual',true,...
    'reversedAxes',3);
  legend off;
  set(h_f,'WindowStyle','Docked');
  if simFluorescence
    h_f.Name = 'Fluorescence light collector illustration';
  else
    h_f.Name = 'Light collector illustration';
  end
  title(h_f.Name);
  box on;grid on;grid minor;

  arrowlength = sqrt((G.nx*G.dx)^2+(G.ny*G.dy)^2+(G.nz*G.dz)^2)/5;
  Zvec = [sin(LC.theta)*cos(LC.phi) , sin(LC.theta)*sin(LC.phi) , cos(LC.theta)];
  Xvec = [sin(LC.phi) , -cos(LC.phi) , 0];
  Yvec = cross(Zvec,Xvec);
  FPC = [LC.x , LC.y , LC.z]; % Focal Plane Center
  FPC_X = FPC + arrowlength*Xvec;
  line([FPC(1) FPC_X(1)],[FPC(2) FPC_X(2)],[FPC(3) FPC_X(3)],'Linewidth',2,'Color','r')
  text(FPC_X(1),FPC_X(2),FPC_X(3),'X','HorizontalAlignment','center','FontSize',18)
  FPC_Y = FPC + arrowlength*Yvec;
  line([FPC(1) FPC_Y(1)],[FPC(2) FPC_Y(2)],[FPC(3) FPC_Y(3)],'Linewidth',2,'Color','r')
  text(FPC_Y(1),FPC_Y(2),FPC_Y(3),'Y','HorizontalAlignment','center','FontSize',18)

  if isfinite(LC.f)
    fieldperimeter = LC.fieldSize/2*(cos(linspace(0,2*pi,100).')*Xvec + sin(linspace(0,2*pi,100).')*Yvec) + FPC;
    h1 = line(fieldperimeter(:,1),fieldperimeter(:,2),fieldperimeter(:,3),'Color','b','LineWidth',2);
    LCC = FPC - Zvec*LC.f; % Light Collector Center
    detectoraperture = LC.diam/2*(cos(linspace(0,2*pi,100).')*Xvec + sin(linspace(0,2*pi,100).')*Yvec) + LCC;
    h2 = line(detectoraperture(:,1),detectoraperture(:,2),detectoraperture(:,3),'Color','r','LineWidth',2);
    legend([h1 h2],'Imaged area','Lens aperture','Location','northeast');
  else
    LCC = FPC;
    detectoraperture = LC.diam/2*(cos(linspace(0,2*pi,100).')*Xvec + sin(linspace(0,2*pi,100).')*Yvec) + LCC;
    h2 = line(detectoraperture(:,1),detectoraperture(:,2),detectoraperture(:,3),'Color','r','LineWidth',2);
    legend(h2,'Fiber aperture','Location','northeast');
  end

  
  if LC.res > 1
    detFraction = 100*mean(mean(sum(MCorFMC.LC.image(:,:,:),3)))*LC.fieldSize^2;
  else
    detFraction = 100*sum(MCorFMC.LC.image(:));
  end
  fprintf(['%.3g%% of ' fluorescenceOrIncident 'light ends up on the detector.\n'],detFraction/P_in);

  if detFraction == 0
    warning('No light was collected on the detector. Are you sure that the light escapes through media that have refractive index 1? Otherwise the photons will not be counted. Depending on your geometry, you could maybe add a layer of air on top of your simulation, or set matchedInterfaces = true, which will set all refractive indices to 1.');
  end

  %% Make NFRdet fluence rate plot
  if ~isscalar(MCorFMC.NFRdet)
    h_f = MCmatlab.NdimSliderPlot(MCorFMC.NFRdet,...
      'nFig',7 + figNumOffset,...
      'axisValues',{G.x,G.y,G.z,MCorFMC.wavelength},...
      'axisLabels',{'x [cm]','y [cm]','z [cm]','λ [nm]',['Normalized fluence rate of collected ' fluorescenceOrIncident 'light [W/cm^2/W.incident]']},...
      'fromZero',true,...
      'axisEqual',true,...
      'reversedAxes',3);
    set(h_f,'WindowStyle','Docked');
    h_f.Name = ['Normalized fluence rate of collected ' fluorescenceOrIncident 'light'];
  end
  
  %% Plot image
  if LC.res > 1 || LC.nTimeBins > 0
    h_f = MCmatlab.NdimSliderPlot(MCorFMC.LC.image,...
      'nFig',8 + figNumOffset,...
      'axisValues',{MCorFMC.LC.X,MCorFMC.LC.Y,MCorFMC.LC.t,MCorFMC.wavelength},...
      'axisLabels',{'X [cm]','Y [cm]','t [s]','λ [nm]','Normalized power [W/W.incident]'},...
      'fromZero',true,...
      'axisDims',[1 2],...
      'axisEqual',true);
    if simFluorescence
      h_f.Name = 'Collected fluorescence light';
    else
      h_f.Name = 'Collected light';
    end
    if LC.res > 1 && LC.nTimeBins > 0
      title({'Normalized time-resolved fluence rate in the image plane','at 1x magnification [W/cm^2/W.incident]'});
      fprintf('Time-resolved light collector data plotted. Note that first time bin includes all\n  photons at earlier times and last time bin includes all photons at later times.\n');
    elseif LC.res > 1
      title({['Normalized ' fluorescenceOrNothing 'fluence rate in the image plane'],' at 1x magnification [W/cm^2/W.incident]'});
    elseif LC.nTimeBins > 0
      title('Normalized time-resolved power on the detector [W/W.incident]');
      fprintf('Time-resolved light collector data plotted. Note that first time bin includes all\n  photons at earlier times and last time bin includes all photons at later times.\n');
    else
      title('Normalized power on the detector [W/W.incident]');
    end
    set(h_f,'WindowStyle','Docked');
  end
end

%% Plot far field of escaped photons
if ~isnan(MCorFMC.farField)
  fprintf(['%.3g%% of ' fluorescenceOrIncident 'light escapes.\n'],100*sum(MCorFMC.farField(:))/P_in);
  if ~isscalar(MCorFMC.farField)
    h_f = plotFarField(9 + figNumOffset,{MCorFMC.wavelength},MCorFMC.farField,...
      {'λ [nm]'});
    if simFluorescence
      h_f.Name = 'Fluorescence far field';
    else
      h_f.Name = 'Far field';
    end
    title(['Far field radiant intensity of escaped ' fluorescenceOrNothing 'photons [W/sr/W.incident]']);
  end
end

%% Plot normalized boundary irradiances
if MCorFMC.boundaryType ~= 0
  plotCuboidSurfaces(10 + figNumOffset,model,simFluorescence);
end

%% Make the NFR plot active so it's the first one people see
if ~isnan(MCorFMC.NFR(1))
  figure(4 + figNumOffset);
end
drawnow;
end

