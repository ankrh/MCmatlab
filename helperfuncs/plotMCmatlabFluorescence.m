function plotMCmatlabFluorescence(FMCinput,FMCoutput)
%   Displays (if calculated)
%       Distribution of fluorescence emitters
%       Absorbed fluorescence power
%       Fluorescence fluence rate of all photons
%       Example paths of some of the photons
%       Distribution of escaped photons in the far field
%       Fluence rates (intensities) on all boundaries
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
h_f = plotVolumetric(12,G.x,G.y,G.z,FluorescenceEmitters,'MCmatlab_fromZero');
h_f.Name = 'Fluorescence emitters';
title('Fluorescence emitter distribution [W/cm^3/W.incident] ')

P_exc_abs = G.dx*G.dy*G.dz*sum(sum(sum(mua_vec(G.M).*FMCinput.MCoutput.F)));
P_flu_emit = G.dx*G.dy*G.dz*sum(sum(sum(FluorescenceEmitters)));
fprintf('\n%.3g%% of absorbed excitation light was re-emitted as fluorescence.\n',100*P_flu_emit/P_exc_abs);

if isfield(FMCoutput,'F')
    %% Make power absorption plot
    mua_vec = [G.mediaProperties_f.mua];
    h_f = plotVolumetric(13,G.x,G.y,G.z,mua_vec(G.M).*FMCoutput.F,'MCmatlab_fromZero');
    h_f.Name = 'Fluorescence power absorption';
    title('Normalized absorbed fluorescence power per unit volume [W/cm^3/W.incident] ')

    %% Make fluence rate plot
    h_f = plotVolumetric(14,G.x,G.y,G.z,FMCoutput.F,'MCmatlab_fromZero');
    h_f.Name = 'Fluorescence fluence rate';
    title('Normalized fluorescence fluence rate (Intensity) [W/cm^2/W.incident] ')

    P_flu_abs = G.dx*G.dy*G.dz*sum(sum(sum(mua_vec(G.M).*FMCoutput.F)));
    fprintf('%.3g%% of emitted fluorescence light was re-absorbed within the cuboid.\n',100*P_flu_abs/P_flu_emit);
end

if isfield(FMCoutput,'ExamplePaths')
    h_f = plotVolumetric(15,G.x,G.y,G.z,G.M,'MCmatlab_GeometryIllustration',G.mediaProperties_f);
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

if(isfield(FMCinput,'LightCollector'))
    %% If there's a fluorescence light collector, show its orientation and the detected light
    h_f = plotVolumetric(16,G.x,G.y,G.z,G.M,'MCmatlab_GeometryIllustration',G.mediaProperties_f);
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
        fprintf('\n%.3g%% of fluorescence light ends up on the detector.\n',100*mean(mean(FMCoutput.Image))*LC.FieldSize^2/P_flu_emit);
    else
        fprintf('\n%.3g%% of fluorescence light ends up on the detector.\n',100*FMCoutput.Image/P_flu_emit);
    end

    if isfield(FMCoutput,'Fdet')
        %% Make Fdet fluence rate plot
        h_f = plotVolumetric(17,G.x,G.y,G.z,FMCoutput.Fdet,'MCmatlab_fromZero');
        h_f.Name = 'Normalized fluence rate of collected fluorescence light';
        title('Normalized fluence rate of collected fluorescence light [W/cm^2/W.incident] ')
    end
    
    if LC.res > 1
        if(~ishandle(18))
            h_f = figure(18);
            h_f.Position = [40 80 1100 650];
        else
            h_f = figure(18);
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

if isfield(FMCoutput,'FarField')
	fprintf('%.3g%% of fluorescence light escapes.\n',100*sum(sum(FMCoutput.FarField))/P_flu_emit);
    farfieldRes = length(FMCoutput.FarField);
    if farfieldRes > 1
        theta_vec = linspace(0,pi,farfieldRes+1).'; % FMCoutput.FFtheta contains the theta values at the centers of the far field pixels, but we will not use those here since we need the corner positions to calculate the solid angle of the pixels
        solidangle_vec = 2*pi*(cos(theta_vec(1:end-1)) - cos(theta_vec(2:end)))/farfieldRes; % Solid angle extended by each far field pixel, as function of theta
        if(~ishandle(19))
            h_f = figure(19);
            h_f.Position = [40 80 1100 650];
        else
            h_f = figure(19);
        end
        h_f.Color = 'w';
		h_f.Name = 'Fluorescence far field';
        clf;
		h_a = axes;

        [ux,uy,minus_uz] = sphere(farfieldRes); % This MATLAB function gives us the correct ux,uy,uz coordinates of the corners of the far field pixels, except that uz has the wrong sign
        uz = -minus_uz;
        surf(ux,uy,uz,FMCoutput.FarField./repmat(solidangle_vec,1,farfieldRes),'EdgeColor','none');

        colormap(inferno);colorbar;
        xlabel('u_x');
        ylabel('u_y');
        zlabel('u_z');
        title('Far field radiant intensity of escaped fluorescence photons [W/sr/W.incident]');
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

if G.boundaryType == 1
    fprintf('%.3g%% of fluorescence light hits the cuboid boundaries.\n',100*(sum(sum((FMCoutput.I_xpos + FMCoutput.I_xneg)*G.dy*G.dz)) + sum(sum((FMCoutput.I_ypos + FMCoutput.I_yneg)*G.dx*G.dz)) + sum(sum((FMCoutput.I_zpos + FMCoutput.I_zneg)*G.dx*G.dy)))/P_flu_emit);
    
    if(~ishandle(20))
        h_f = figure(20);
        h_f.Position = [40 80 1100 650];
    else
        h_f = figure(20);
    end
    h_f.Color = 'w';
    h_f.Name = 'Boundary fluorescence fluence rate';
    clf;
    h_a = axes;
    title('Boundary fluorescence fluence rate [W/cm^2/W.incident]');
    
    x = round([(G.x - G.dx/2) , (max(G.x) + G.dx/2)],15);
    y = round([(G.y - G.dy/2) , (max(G.y) + G.dy/2)],15);
    z = round([(G.z - G.dz/2) , (max(G.z) + G.dz/2)],15);
    xl = x(1); % x low
    xh = x(end); % x high
    yl = y(1); % y low
    yh = y(end); % y high
    zl = z(1); % z low
    zh = z(end); % z high
    I_xneg_pad = FMCoutput.I_xneg;
    I_xneg_pad(G.ny+1,G.nz+1) = 0;
    I_yneg_pad = FMCoutput.I_yneg;
    I_yneg_pad(G.nx+1,G.nz+1) = 0;
    I_zneg_pad = FMCoutput.I_zneg;
    I_zneg_pad(G.nx+1,G.ny+1) = 0;
    I_xpos_pad = FMCoutput.I_xpos;
    I_xpos_pad(G.ny+1,G.nz+1) = 0;
    I_ypos_pad = FMCoutput.I_ypos;
    I_ypos_pad(G.nx+1,G.nz+1) = 0;
    I_zpos_pad = FMCoutput.I_zpos;
    I_zpos_pad(G.nx+1,G.ny+1) = 0;

    surface(repmat(xl,G.ny+1,G.nz+1),repmat(y',     1,G.nz+1),repmat(z ,G.ny+1,     1),I_xneg_pad,'LineStyle','none');
    surface(repmat(x',     1,G.nz+1),repmat(yl,G.nx+1,G.nz+1),repmat(z ,G.nx+1,     1),I_yneg_pad,'LineStyle','none');
    surface(repmat(x',     1,G.ny+1),repmat(y ,G.nx+1,     1),repmat(zl,G.nx+1,G.ny+1),I_zneg_pad,'LineStyle','none');
    surface(repmat(xh,G.ny+1,G.nz+1),repmat(y',     1,G.nz+1),repmat(z ,G.ny+1,     1),I_xpos_pad,'LineStyle','none');
    surface(repmat(x',     1,G.nz+1),repmat(yh,G.nx+1,G.nz+1),repmat(z ,G.nx+1,     1),I_ypos_pad,'LineStyle','none');
    surface(repmat(x',     1,G.ny+1),repmat(y ,G.nx+1,     1),repmat(zh,G.nx+1,G.ny+1),I_zpos_pad,'LineStyle','none');
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
elseif G.boundaryType == 2
    fprintf('%.3g%% of fluorescence light hits the top cuboid boundary.\n',100*(sum(sum(FMCoutput.I_zneg*G.dx*G.dy)))/P_flu_emit);
    
    if(~ishandle(20))
        h_f = figure(20);
        h_f.Position = [40 80 1100 650];
    else
        h_f = figure(20);
    end
    h_f.Color = 'w';
    h_f.Name = 'Boundary fluorescence fluence rate';
    clf;
    imagesc(size(FMCoutput.I_zneg,1)/2*[-G.dx G.dx],size(FMCoutput.I_zneg,2)/2*[-G.dy G.dy],FMCoutput.I_zneg.');
	line(G.Lx/2*[-1 -1 1 1 -1],G.Ly/2*[-1 1 1 -1 -1],'Color',[1 1 1],'Linestyle','--');
    set(gca,'YDir','normal');
    title('Boundary fluorescence fluence rate [W/cm^2/W.incident]');
    colormap(inferno);
    colorbar;
    axis tight
    axis equal
    xlabel('x [cm]');
    ylabel('y [cm]');
    set(gca,'fontsize',18)
end
drawnow;
end
