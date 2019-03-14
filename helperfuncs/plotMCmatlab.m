function plotMCmatlab(MCinput,MCoutput)
%   Displays (if calculated)
%       Absorbed power
%       Fluence rate of all photons
%       Example paths of some of the photons
%       Distribution of escaped photons in the far field
%       Fluence rates (intensities) on all boundaries
%	And, if a light collector was defined, displays (if calculated)
%		An illustration of the light collector angle and placement
%		Image generated (which might be time-resolved)
%       Fluence rate of photons that hit the light collector
%
%   Requires
%       plotVolumetric.m
%
%	See also runMonteCarlo

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

G = MCinput.G;
    
if isfield(MCoutput,'F')
    %% Make power absorption plot
    mua_vec = [G.mediaProperties.mua];
    h_f = plotVolumetric(4,G.x,G.y,G.z,mua_vec(G.M).*MCoutput.F,'MCmatlab_fromZero');
    h_f.Name = 'Normalized power absorption';
    title('Normalized absorbed power per unit volume [W/cm^3/W.incident] ')

    %% Make fluence rate plot
    h_f = plotVolumetric(5,G.x,G.y,G.z,MCoutput.F,'MCmatlab_fromZero');
    h_f.Name = 'Normalized fluence rate';
    title('Normalized fluence rate (Intensity) [W/cm^2/W.incident] ')

    fprintf('\n%.3g%% of incident light was absorbed within the cuboid.\n',100*G.dx*G.dy*G.dz*sum(sum(sum(mua_vec(G.M).*MCoutput.F))));
end

if isfield(MCoutput,'ExamplePaths')
    h_f = plotVolumetric(6,G.x,G.y,G.z,G.M,'MCmatlab_GeometryIllustration',G.mediaProperties);
    h_f.Name = 'Photon paths';
    title('Photon paths');
    box on;grid on;grid minor;

    previousNaNidx = 1;
    linenumber = 0;
    for idx=3:size(MCoutput.ExamplePaths,2)
        if isnan(MCoutput.ExamplePaths(1,idx))
            linenumber = linenumber + 1;
            xdata = MCoutput.ExamplePaths(1,previousNaNidx+1:idx-1);
            ydata = MCoutput.ExamplePaths(2,previousNaNidx+1:idx-1);
            zdata = MCoutput.ExamplePaths(3,previousNaNidx+1:idx-1);
            adata = MCoutput.ExamplePaths(4,previousNaNidx+1:idx-1);
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

if(isfield(MCinput,'LightCollector'))
    %% If there's a light collector, show its orientation and the detected light
    h_f = plotVolumetric(7,G.x,G.y,G.z,G.M,'MCmatlab_GeometryIllustration',G.mediaProperties);
    h_f.Name = 'Light collector illustration';
    title('Light collector illustration');
    box on;grid on;grid minor;

    LC = MCinput.LightCollector;
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
        fprintf('\n%.3g%% of incident light ends up on the detector.\n',100*mean(mean(sum(MCoutput.Image,3)))*LC.FieldSize^2);
    else
        fprintf('\n%.3g%% of incident light ends up on the detector.\n',100*sum(MCoutput.Image,3));
    end

    if isfield(MCoutput,'Fdet')
        %% Make Fdet fluence rate plot
        h_f = plotVolumetric(8,G.x,G.y,G.z,MCoutput.Fdet,'MCmatlab_fromZero');
        h_f.Name = 'Normalized fluence rate of collected light';
        title('Normalized fluence rate of collected light [W/cm^2/W.incident] ')
    end
    
    if(~isfield(LC,'nTimeBins'))
        LC.nTimeBins = 0;
    end
    
    if LC.nTimeBins > 0
        timevector = (-1/2:(LC.nTimeBins+1/2))*(LC.tEnd-LC.tStart)/LC.nTimeBins + LC.tStart;
    end

    if LC.res > 1 && LC.nTimeBins > 0
        h_f = plotVolumetric(9,MCoutput.X,MCoutput.Y,timevector,MCoutput.Image,'slicePositions',[1 1 0]);
        h_f.Name = 'Image';
        xlabel('X [cm]');
        ylabel('Y [cm]');
        zlabel('Time [s]');
        title({'Normalized time-resolved fluence rate in the image plane','at 1x magnification [W/cm^2/W.incident]'});
        fprintf('Time-resolved light collector data plotted. Note that first time bin includes all\n  photons at earlier times and last time bin includes all photons at later times.\n');
    elseif LC.res > 1
        if(~ishandle(9))
            h_f = figure(9);
            h_f.Position = [40 80 1100 650];
        else
            h_f = figure(9);
        end
		h_f.Color = 'w';
        clf;
        h_f.Name = 'Image';
        imagesc(MCoutput.X,MCoutput.Y,MCoutput.Image.');
        title({'Normalized fluence rate in the image plane',' at 1x magnification [W/cm^2/W.incident]'});
        axis xy;axis equal;axis tight;xlabel('X [cm]');ylabel('Y [cm]');
        set(gca,'FontSize',18);
        colormap(inferno);
        colorbar;
    elseif LC.nTimeBins > 0
        if(~ishandle(9))
            h_f = figure(9);
            h_f.Position = [40 80 1100 650];
        else
            h_f = figure(9);
        end
		h_f.Color = 'w';
        clf;
        h_b = bar(timevector,squeeze(MCoutput.Image),1,'FaceColor','flat');
        h_b.CData(1  ,:) = [.5 0 .5];
        h_b.CData(end,:) = [.5 0 .5];
        title('Normalized time-resolved power on the detector');
        xlabel('Time [s]'); ylabel('Normalized power [W/W.incident]'); grid on; grid minor;
        set(gca,'FontSize',18);
        fprintf('Time-resolved light collector data plotted. Note that first time bin includes all\n  photons at earlier times and last time bin includes all photons at later times.\n');
    end
end

if isfield(MCoutput,'FarField')
	fprintf('%.3g%% of incident light escapes.\n',100*sum(sum(MCoutput.FarField)));
    farfieldRes = length(MCoutput.FarField);
    if farfieldRes > 1
        theta_vec = linspace(0,pi,farfieldRes+1).'; % MCoutput.FFtheta contains the theta values at the centers of the far field pixels, but we will not use those here since we need the corner positions to calculate the solid angle of the pixels
        solidangle_vec = 2*pi*(cos(theta_vec(1:end-1)) - cos(theta_vec(2:end)))/farfieldRes; % Solid angle extended by each far field pixel, as function of theta
        if(~ishandle(10))
            h_f = figure(10);
            h_f.Position = [40 80 1100 650];
        else
            h_f = figure(10);
        end
        h_f.Color = 'w';
		h_f.Name = 'Far field';
        clf;
		h_a = axes;

        [ux,uy,minus_uz] = sphere(farfieldRes); % This MATLAB function gives us the correct ux,uy,uz coordinates of the corners of the far field pixels, except that uz has the wrong sign
        uz = -minus_uz;
        surf(ux,uy,uz,MCoutput.FarField./repmat(solidangle_vec,1,farfieldRes),'EdgeColor','none');

        colormap(inferno);colorbar;
        xlabel('u_x');
        ylabel('u_y');
        zlabel('u_z');
        title('Far field radiant intensity of escaped photons [W/sr/W.incident]');
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
    fprintf('%.3g%% of incident light hits the cuboid boundaries.\n',100*(sum(sum((MCoutput.I_xpos + MCoutput.I_xneg)*G.dy*G.dz)) + sum(sum((MCoutput.I_ypos + MCoutput.I_yneg)*G.dx*G.dz)) + sum(sum((MCoutput.I_zpos + MCoutput.I_zneg)*G.dx*G.dy))));
    
    if(~ishandle(11))
        h_f = figure(11);
        h_f.Position = [40 80 1100 650];
    else
        h_f = figure(11);
    end
    h_f.Color = 'w';
    h_f.Name = 'Boundary fluence rate';
    clf;
    h_a = axes;
    title('Boundary fluence rate [W/cm^2/W.incident]');
    
    x = round([(G.x - G.dx/2) , (max(G.x) + G.dx/2)],15);
    y = round([(G.y - G.dy/2) , (max(G.y) + G.dy/2)],15);
    z = round([(G.z - G.dz/2) , (max(G.z) + G.dz/2)],15);
    xl = x(1); % x low
    xh = x(end); % x high
    yl = y(1); % y low
    yh = y(end); % y high
    zl = z(1); % z low
    zh = z(end); % z high
    I_xneg_pad = MCoutput.I_xneg;
    I_xneg_pad(G.ny+1,G.nz+1) = 0;
    I_yneg_pad = MCoutput.I_yneg;
    I_yneg_pad(G.nx+1,G.nz+1) = 0;
    I_zneg_pad = MCoutput.I_zneg;
    I_zneg_pad(G.nx+1,G.ny+1) = 0;
    I_xpos_pad = MCoutput.I_xpos;
    I_xpos_pad(G.ny+1,G.nz+1) = 0;
    I_ypos_pad = MCoutput.I_ypos;
    I_ypos_pad(G.nx+1,G.nz+1) = 0;
    I_zpos_pad = MCoutput.I_zpos;
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
    if ~(isfield(MCinput.Beam,'beamType') && MCinput.Beam.beamType == 2) 
        fprintf('%.3g%% of incident light hits the top cuboid boundary.\n',100*(sum(sum(MCoutput.I_zneg*G.dx*G.dy))));
    end
    
    if(~ishandle(11))
        h_f = figure(11);
        h_f.Position = [40 80 1100 650];
    else
        h_f = figure(11);
    end
    h_f.Color = 'w';
    h_f.Name = 'Boundary fluence rate';
    clf;
    imagesc(size(MCoutput.I_zneg,1)/2*[-G.dx G.dx],size(MCoutput.I_zneg,2)/2*[-G.dy G.dy],MCoutput.I_zneg.');
	line(G.Lx/2*[-1 -1 1 1 -1],G.Ly/2*[-1 1 1 -1 -1],'Color',[1 1 1],'Linestyle','--');
    set(gca,'YDir','normal');
    title('Boundary fluence rate [W/cm^2/W.incident]');
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

