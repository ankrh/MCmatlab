function plotMCmatlabMC(model,varargin)
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
%       plotVolumetric.m
%
%   See also runMonteCarlo

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

% 1: Geometry illustration
% 2: Media properties
% 3: Absorption
% 4: NFR
% 5: Paths
% 6: LC illustration
% 7: NFRdet
% 8: Image
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

LC = MCorFMC.LC;


%% Plot emitter distribution
P_in = 1;
if simFluorescence
  h_f = plotVolumetric.plotVolumetric(11,G.x,G.y,G.z,model.FMC.sourceDistribution,'MCmatlab_fromZero');
  set(h_f,'WindowStyle','Docked');
  h_f.Name = 'Fluorescence emitters';
  title('Fluorescence emitter distribution [W/cm^3/W.incident]')

  % Calculate 3D absorption distribution, which may be FR or T dependent
  mua_3d = NaN(size(G.M_raw));
  for i=1:length(model.MC.mediaProperties_funcHandles)
    if isnan(model.MC.FR(1))
      FR = zeros(size(G.M_raw));
    else
      FR = model.MC.FR;
    end
    if isnan(model.HS.T(1))
      T = zeros(size(G.M_raw));
    else
      T = model.HS.T;
    end
    if isa(model.MC.mediaProperties_funcHandles(i).mua,'function_handle')
      mua_3d(G.M_raw == i) = model.MC.mediaProperties_funcHandles(i).mua(FR(G.M_raw == i),T(G.M_raw == i));
    else
      mua_3d(G.M_raw == i) = model.MC.mediaProperties_funcHandles(i).mua;
    end
  end
  P_exc_abs = G.dx*G.dy*G.dz*sum(sum(sum(mua_3d.*model.MC.NFR)));
  clear mua_3d T
  P_flu_emit = G.dx*G.dy*G.dz*sum(sum(sum(model.FMC.sourceDistribution)));
  fprintf('%.3g%% of absorbed incident light was re-emitted as fluorescence.\n',100*P_flu_emit/P_exc_abs);
  P_in = P_flu_emit;
end

%% Make media properties plot
h_f = plotMediaProperties(2+figNumOffset,model,simType);
set(h_f,'WindowStyle','Docked');
if simFluorescence
  h_f.Name = 'Fluorescence media properties';
else
  h_f.Name = 'Media properties';
end

if ~isnan(MCorFMC.NFR(1))
  %% Make power absorption plot
  % Calculate 3D absorption distribution, which may be FR or T dependent
  mua_vec = [MCorFMC.mediaProperties.mua];
  h_f = plotVolumetric.plotVolumetric(3 + figNumOffset,G.x,G.y,G.z,mua_vec(MCorFMC.M).*MCorFMC.NFR,'MCmatlab_fromZero');
  set(h_f,'WindowStyle','Docked');
  h_f.Name = ['Normalized ' fluorescenceOrNothing 'absorption'];
  title(['Normalized ' fluorescenceOrNothing 'absorbed power per unit volume [W/cm^3/W.incident]'])
  
  %% Make fluence rate plot
  h_f = plotVolumetric.plotVolumetric(4 + figNumOffset,G.x,G.y,G.z,MCorFMC.NFR,'MCmatlab_fromZero');
  set(h_f,'WindowStyle','Docked');
  h_f.Name = ['Normalized ' fluorescenceOrNothing 'fluence rate'];
  title(['Normalized ' fluorescenceOrNothing 'fluence rate [W/cm^2/W.incident]'])
  
  fprintf(['%.3g%% of ' fluorescenceOrIncident 'light was absorbed within the cuboid.\n'],100*G.dx*G.dy*G.dz*sum(sum(sum(mua_vec(MCorFMC.M).*MCorFMC.NFR)))/P_in);
end

%% Plot example paths
if MCorFMC.nExamplePaths > 0
  h_f = plotVolumetric.plotVolumetric(5 + figNumOffset,G.x,G.y,G.z,G.M_raw,'MCmatlab_GeometryIllustration',MCorFMC.mediaProperties_funcHandles);
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
              'CDataMapping','direct',...
              'AlphaDataMapping','direct',...
              'LineWidth',2);
    end
  end
end

if MCorFMC.useLightCollector
  %% If there's a light collector, show its orientation and the detected light
  h_f = plotVolumetric.plotVolumetric(6 + figNumOffset,G.x,G.y,G.z,G.M_raw,'MCmatlab_GeometryIllustration',MCorFMC.mediaProperties_funcHandles);
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
    detFraction = 100*mean(mean(sum(MCorFMC.LC.image,3)))*LC.fieldSize^2;
  else
    detFraction = 100*sum(MCorFMC.LC.image,3);
  end
  fprintf(['%.3g%% of ' fluorescenceOrIncident 'light ends up on the detector.\n'],detFraction/P_in);

  if detFraction == 0
    warning('No light was collected on the detector. Are you sure that the light escapes through media that have refractive index 1? Otherwise the photons will not be counted. Depending on your geometry, you could maybe add a layer of air on top of your simulation, or set matchedInterfaces = true, which will set all refractive indices to 1.');
  end

  if ~isnan(MCorFMC.NFRdet(1))
    %% Make NFRdet fluence rate plot
    h_f = plotVolumetric.plotVolumetric(7 + figNumOffset,G.x,G.y,G.z,MCorFMC.NFRdet,'MCmatlab_fromZero');
    set(h_f,'WindowStyle','Docked');
    h_f.Name = ['Normalized fluence rate of collected ' fluorescenceOrIncident 'light'];
    title(['Normalized fluence rate of collected ' fluorescenceOrIncident 'light [W/cm^2/W.incident]'])
  end
  
  if ~simFluorescence && LC.nTimeBins > 0
    timevector = (-1/2:(LC.nTimeBins+1/2))*(LC.tEnd-LC.tStart)/LC.nTimeBins + LC.tStart;
  end

  if LC.res > 1 && ~simFluorescence && LC.nTimeBins > 0
    h_f = plotVolumetric.plotVolumetric(8 + figNumOffset,LC.X,LC.Y,timevector,LC.image,'slicePositions',[1 1 0]);
    h_f.Name = 'Image';
    xlabel('X [cm]');
    ylabel('Y [cm]');
    zlabel('Time [s]');
    title({'Normalized time-resolved fluence rate in the image plane','at 1x magnification [W/cm^2/W.incident]'});
    fprintf('Time-resolved light collector data plotted. Note that first time bin includes all\n  photons at earlier times and last time bin includes all photons at later times.\n');
  elseif LC.res > 1
    h_f = figure(8 + figNumOffset);
    h_f.Color = 'w';
    clf;
    if simFluorescence
      h_f.Name = 'Fluorescence image';
    else
      h_f.Name = 'Image';
    end
    imagesc(LC.X,LC.Y,LC.image.');
    title({['Normalized ' fluorescenceOrNothing 'fluence rate in the image plane'],' at 1x magnification [W/cm^2/W.incident]'});
    axis xy;axis equal;axis tight;xlabel('X [cm]');ylabel('Y [cm]');
    set(gca,'FontSize',18);
    colormap(inferno);
    colorbar;
  elseif ~simFluorescence && LC.nTimeBins > 0
    h_f = figure(8 + figNumOffset);
    h_f.Color = 'w';
    clf;
    h_f.Name = 'Time-resolved detected power';
    h_b = bar(timevector,squeeze(LC.image),1,'FaceColor','flat');
    h_b.CData(1  ,:) = [.5 0 .5];
    h_b.CData(end,:) = [.5 0 .5];
    title('Normalized time-resolved power on the detector');
    xlabel('Time [s]'); ylabel('Normalized power [W/W.incident]'); grid on; grid minor;
    set(gca,'FontSize',18);
    fprintf('Time-resolved light collector data plotted. Note that first time bin includes all\n  photons at earlier times and last time bin includes all photons at later times.\n');
  end
  set(h_f,'WindowStyle','Docked');
end

if ~isnan(MCorFMC.farField)
  fprintf(['%.3g%% of ' fluorescenceOrIncident 'light escapes.\n'],100*sum(sum(MCorFMC.farField))/P_in);
  farFieldRes = length(MCorFMC.farField);
  if farFieldRes > 1
    theta_vec = linspace(0,pi,farFieldRes+1).'; % MCorFMC.farFieldTheta contains the theta values at the centers of the far field pixels, but we will not use those here since we need the corner positions to calculate the solid angle of the pixels
    solidAngle_vec = 2*pi*(cos(theta_vec(1:end-1)) - cos(theta_vec(2:end)))/farFieldRes; % Solid angle extended by each far field pixel, as function of theta
    h_f = figure(9 + figNumOffset);
    set(h_f,'WindowStyle','Docked');
    h_f.Color = 'w';
    if simFluorescence
      h_f.Name = 'Fluorescence far field';
    else
      h_f.Name = 'Far field';
    end
    clf;
    h_a = axes;

    [ux,uy,minus_uz] = sphere(farFieldRes); % This MATLAB function gives us the correct ux,uy,uz coordinates of the corners of the far field pixels, except that uz has the wrong sign
    uz = -minus_uz;
    plotdata = MCorFMC.farField./repmat(solidAngle_vec,1,farFieldRes);
    plotdata(1,:) = mean(plotdata(1,:));
    plotdata(end,:) = mean(plotdata(end,:));
    surf(ux,uy,uz,plotdata,'EdgeColor','none');

    colormap(inferno);colorbar;
    xlabel('u_x');
    ylabel('u_y');
    zlabel('u_z');
    title(['Far field radiant intensity of escaped ' fluorescenceOrNothing 'photons [W/sr/W.incident]']);
    set(gca,'FontSize',18);
    axis equal;
    xlim([-1 1]);
    ylim([-1 1]);
    zlim([-1 1]);
    h_a.ZDir = 'reverse';
    h_a.CLim = [0 h_a.CLim(2)];
    rotate3d on
    if ~verLessThan('matlab','9.0')
      setAxes3DPanAndZoomStyle(zoom(gca),gca,'camera');
    end
  end
end

%% Plot normalized boundary irradiances
if MCorFMC.boundaryType == 1
  fprintf(['%.3g%% of ' fluorescenceOrIncident 'light hits the cuboid boundaries.\n'],100*(sum(sum((MCorFMC.NI_xpos + MCorFMC.NI_xneg)*G.dy*G.dz)) + sum(sum((MCorFMC.NI_ypos + MCorFMC.NI_yneg)*G.dx*G.dz)) + sum(sum((MCorFMC.NI_zpos + MCorFMC.NI_zneg)*G.dx*G.dy)))/P_in);
  
  h_f = figure(10 + figNumOffset);
  set(h_f,'WindowStyle','Docked');
  h_f.Color = 'w';
  h_f.Name = ['Normalized ' fluorescenceOrNothing 'boundary irradiance'];
  clf;
  title(['Normalized ' fluorescenceOrNothing 'boundary irradiance [W/cm^2/W.incident]']);
  
  x = round([(G.x - G.dx/2) , (max(G.x) + G.dx/2)],15);
  y = round([(G.y - G.dy/2) , (max(G.y) + G.dy/2)],15);
  z = round([(G.z - G.dz/2) , (max(G.z) + G.dz/2)],15);
  xl = x(1); % x low
  xh = x(end); % x high
  yl = y(1); % y low
  yh = y(end); % y high
  zl = z(1); % z low
  zh = z(end); % z high
  NI_xneg_pad = MCorFMC.NI_xneg;
  NI_xneg_pad(G.ny+1,G.nz+1) = 0;
  NI_yneg_pad = MCorFMC.NI_yneg;
  NI_yneg_pad(G.nx+1,G.nz+1) = 0;
  NI_zneg_pad = MCorFMC.NI_zneg;
  NI_zneg_pad(G.nx+1,G.ny+1) = 0;
  NI_xpos_pad = MCorFMC.NI_xpos;
  NI_xpos_pad(G.ny+1,G.nz+1) = 0;
  NI_ypos_pad = MCorFMC.NI_ypos;
  NI_ypos_pad(G.nx+1,G.nz+1) = 0;
  NI_zpos_pad = MCorFMC.NI_zpos;
  NI_zpos_pad(G.nx+1,G.ny+1) = 0;

  surface(repmat(xl,G.ny+1,G.nz+1),repmat(y',     1,G.nz+1),repmat(z ,G.ny+1,     1),NI_xneg_pad,'LineStyle','none');
  surface(repmat(x',     1,G.nz+1),repmat(yl,G.nx+1,G.nz+1),repmat(z ,G.nx+1,     1),NI_yneg_pad,'LineStyle','none');
  surface(repmat(x',     1,G.ny+1),repmat(y ,G.nx+1,     1),repmat(zl,G.nx+1,G.ny+1),NI_zneg_pad,'LineStyle','none');
  surface(repmat(xh,G.ny+1,G.nz+1),repmat(y',     1,G.nz+1),repmat(z ,G.ny+1,     1),NI_xpos_pad,'LineStyle','none');
  surface(repmat(x',     1,G.nz+1),repmat(yh,G.nx+1,G.nz+1),repmat(z ,G.nx+1,     1),NI_ypos_pad,'LineStyle','none');
  surface(repmat(x',     1,G.ny+1),repmat(y ,G.nx+1,     1),repmat(zh,G.nx+1,G.ny+1),NI_zpos_pad,'LineStyle','none');
  line(G.Lx/2*[1 1 1 1 -1 -1 -1 -1 1],G.Ly/2*[1 1 -1 -1 -1 -1 1 1 1],G.Lz*[1 0 0 1 1 0 0 1 1],'Color',[0.5 0.5 0.5]);
  line(G.Lx/2*[1 1 1 1 -1 -1 -1 -1 1],G.Ly/2*[1 1 -1 -1 -1 -1 1 1 1],G.Lz*[0 1 1 0 0 1 1 0 0],'Color',[0.5 0.5 0.5]);

  set(gca,'ZDir','reverse');
  colormap(inferno);
  colorbar;
  axis tight
  axis equal
  xlabel('x [cm]');
  ylabel('y [cm]');
  zlabel('z [cm]');
  set(gca,'fontsize',18)
  view(3)
  rotate3d on
  if ~verLessThan('matlab','9.0')
    setAxes3DPanAndZoomStyle(zoom(gca),gca,'camera');
  end
elseif MCorFMC.boundaryType == 2
  if simFluorescence || isnan(MCorFMC.lightSource.sourceType) || MCorFMC.lightSource.sourceType ~= 2
    fprintf(['%.3g%% of ' fluorescenceOrIncident 'light hits the top cuboid boundary.\n'],100*(sum(sum(MCorFMC.NI_zneg*G.dx*G.dy)))/P_in);
  end
  
  h_f = figure(10 + figNumOffset);
  set(h_f,'WindowStyle','Docked');
  h_f.Color = 'w';
  h_f.Name = ['Normalized ' fluorescenceOrNothing 'boundary irradiance'];
  clf;
  imagesc(size(MCorFMC.NI_zneg,1)/2*[-G.dx G.dx],size(MCorFMC.NI_zneg,2)/2*[-G.dy G.dy],MCorFMC.NI_zneg.');
  line(G.Lx/2*[-1 -1 1 1 -1],G.Ly/2*[-1 1 1 -1 -1],'Color',[1 1 1],'Linestyle','--');
  set(gca,'YDir','normal');
  title(['Normalized ' fluorescenceOrNothing 'boundary irradiance [W/cm^2/W.incident]']);
  colormap(inferno);
  colorbar;
  axis tight
  axis equal
  xlabel('x [cm]');
  ylabel('y [cm]');
  set(gca,'fontsize',18)
end
if ~isnan(MCorFMC.NFR(1))
  figure(4 + figNumOffset); % Make the NFR plot active so it's the first one people see
end
drawnow;
end

