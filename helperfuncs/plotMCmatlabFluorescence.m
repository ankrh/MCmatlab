function plotMCmatlabFluorescence(FMCinput,FMCoutput)
%   Displays (if calculated)
%       Distribution of fluorescence emitters
%       Absorbed fluorescence power
%       Fluorescence fluence rate of all photons
%       Example paths of some of the photons
%	And, if a light collector was defined, displays (if calculated)
%		An illustration of the light collector angle and placement
%		Image generated
%       Fluorescence fluence rate of photons that hit the light collector
%
%   Requires
%       plotVolumetric.m
%
%	See also runMonteCarloFluorescence

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

G = FMCinput.G;
mua_vec = [G.mediaProperties.mua];

%% Plot emitter distribution
Y_vec = [G.mediaProperties.Y]; % The media's fluorescence power efficiencies
FluorescenceEmitters = Y_vec(G.M).*mua_vec(G.M).*FMCinput.MCoutput.F; % [W/cm^3]
h_f = plotVolumetric(10,G.x,G.y,G.z,FluorescenceEmitters,'MCmatlab_fromZero');
h_f.Name = 'Fluorescence emitters';
title('Fluorescence emitter distribution [W/cm^3/W.incident] ')

P_exc_abs = G.dx*G.dy*G.dz*sum(sum(sum(mua_vec(G.M).*FMCinput.MCoutput.F)));
P_flu_emit = G.dx*G.dy*G.dz*sum(sum(sum(FluorescenceEmitters)));
fprintf('\n%.3g%% of absorbed excitation light was re-emitted as fluorescence.\n',100*P_flu_emit/P_exc_abs);

if isfield(FMCoutput,'F')
    %% Make power absorption plot
    mua_vec = [G.mediaProperties_f.mua];
    h_f = plotVolumetric(11,G.x,G.y,G.z,mua_vec(G.M).*FMCoutput.F,'MCmatlab_fromZero');
    h_f.Name = 'Fluorescence power absorption';
    title('Normalized absorbed fluorescence power per unit volume [W/cm^3/W.incident] ')

    %% Make fluence rate plot
    h_f = plotVolumetric(12,G.x,G.y,G.z,FMCoutput.F,'MCmatlab_fromZero');
    h_f.Name = 'Fluorescence fluence rate';
    title('Normalized fluorescence fluence rate (Intensity) [W/cm^2/W.incident] ')

    P_flu_abs = G.dx*G.dy*G.dz*sum(sum(sum(mua_vec(G.M).*FMCoutput.F)));
    fprintf('%.3g%% of emitted fluorescence light was re-absorbed within the cuboid.\n',100*P_flu_abs/P_flu_emit);
end

if isfield(FMCoutput,'ExamplePaths')
    h_f = plotVolumetric(13,G.x,G.y,G.z,G.M,'MCmatlab_GeometryIllustration',G.mediaProperties_f);
    h_f.Name = 'Fluorescence photon paths';
    title('Fluorescence photon paths');
    box on;grid on;grid minor;

    previousNaNidx = 1;
    linenumber = 0;
    for idx=3:size(FMCoutput.ExamplePaths,2)
        if isnan(FMCoutput.ExamplePaths(1,idx))
            linenumber = linenumber + 1;
            xdata = FMCoutput.ExamplePaths(1,previousNaNidx+1:idx-1);
            ydata = FMCoutput.ExamplePaths(2,previousNaNidx+1:idx-1);
            zdata = FMCoutput.ExamplePaths(3,previousNaNidx+1:idx-1);
            adata = FMCoutput.ExamplePaths(4,previousNaNidx+1:idx-1);
            previousNaNidx = idx;
            surface([xdata;xdata],...
				    [ydata;ydata],...
				    [zdata;zdata],...
                    'AlphaData',[adata;adata],...
				    'EdgeColor',[0 0 0],...
                    'EdgeAlpha','flat',...
                    'FaceColor','none',...
                    'CDataMapping','direct',...
                    'LineWidth',2);
        end
    end
end

if(isfield(FMCinput,'LightCollector'))
    %% If there's a fluorescence light collector, show its orientation and the detected light
    h_f = plotVolumetric(14,G.x,G.y,G.z,G.M,'MCmatlab_GeometryIllustration',G.mediaProperties_f);
    h_f.Name = 'Fluorescence light collector illustration';
    title('Fluorescence light collector illustration');
    box on;grid on;grid minor;

    LC = FMCinput.LightCollector;
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
        fieldperimeter = LC.FieldSize/2*(cos(linspace(0,2*pi,100).')*Xvec + sin(linspace(0,2*pi,100).')*Yvec) + FPC;
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
        fprintf('\n%.3g%% of fluorescence light ends up on the detector.\n',100*mean(mean(FMCoutput.Image))*LC.FieldSize^2);
    else
        fprintf('\n%.3g%% of fluorescence light ends up on the detector.\n',100*FMCoutput.Image);
    end

    if isfield(FMCoutput,'Fdet')
        %% Make Fdet fluence rate plot
        h_f = plotVolumetric(15,G.x,G.y,G.z,FMCoutput.Fdet,'MCmatlab_fromZero');
        h_f.Name = 'Normalized fluence rate of collected fluorescence light';
        title('Normalized fluence rate of collected fluorescence light [W/cm^2/W.incident] ')
    end
    
    if LC.res > 1
        if(~ishandle(16))
            h_f = figure(16);
            h_f.Position = [40 80 1100 650];
        else
            h_f = figure(16);
        end
		h_f.Color = 'w';
        clf;
        h_f.Name = 'Fluorescence image';
        imagesc(FMCoutput.X,FMCoutput.Y,FMCoutput.Image.');
        title({'Normalized fluence rate in the fluorescence image plane',' at 1x magnification [W/cm^2/W.incident]'});axis xy;axis equal;axis tight;xlabel('X [cm]');ylabel('Y [cm]');
        set(gca,'FontSize',18);
        colormap(inferno);
        colorbar;
    end
end
end
