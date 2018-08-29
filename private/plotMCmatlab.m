function plotMCmatlab(MCinput,MCoutput)
%   Created 2018 by Dominik Marti and Anders K. Hansen, DTU Fotonik
%
%   This function was inspired by lookmcxyz.m of the mcxyz MC program hosted at omlc.org
%
%   Input
%       name
%           the basename of the data files
%
%   Displays
%       Geometry cuboid
%       Media optical, thermal and fluorescence properties
%   If Monte Carlo output data exists, displays
%       Absorbed power
%       Fluence rate
%   And, if fluorescence Monte Carlo output data exists, displays
%       Media optical and thermal properties at the fluorescence wavelength
%       Distribution of fluorescence emitters
%       Absorbed fluorescence power
%       Fluorescence fluence rate
%	And, if light collectors were defined, displays
%		An illustration of the light collector angle and placement
%		Image generated
%
%   Requires
%       plotVolumetric.m
%       plotMediaProperties.m
%

G = MCinput.G;

%% Make power absorption plot
mua_vec = [G.mediaProperties.mua];
h_f = plotVolumetric(4,G.x,G.y,G.z,mua_vec(G.M).*MCoutput.F,'MCmatlab_fromZero');
h_f.Name = 'Normalized power absorption';
title('Normalized absorbed power per unit volume [W/cm^3/W.incident] ')

%% Make fluence rate plot
h_f = plotVolumetric(5,G.x,G.y,G.z,MCoutput.F,'MCmatlab_fromZero');
h_f.Name = 'Normalized fluence rate';
title('Normalized fluence rate (Intensity) [W/cm^2/W.incident] ')

fprintf('\n%.3g%% of the input light was absorbed within the cuboid.\n',100*G.dx*G.dy*G.dz*sum(sum(sum(mua_vec(G.M).*MCoutput.F))));

if(MCinput.useLightCollector)
    %% If there's a light collector, show its orientation and the detected light
    h_f = plotVolumetric(6,G.x,G.y,G.z,G.M,'MCmatlab_GeometryIllustration',G.mediaProperties);
    h_f.Name = 'Light collector illustration';
    title('Light collector illustration');
    box on;grid on;grid minor;

    LC = MCinput.LightCollector;
    arrowlength = sqrt((G.nx*G.dx)^2+(G.ny*G.dy)^2+(G.nz*G.dz)^2)/5;
    Zvec = [sin(LC.theta)*cos(LC.phi) , sin(LC.theta)*sin(LC.phi) , cos(LC.theta)];
    Xvec = [sin(LC.phi) , -cos(LC.phi) , 0];
    Yvec = cross(Zvec,Xvec);
    FPC = [LC.xFPC , LC.yFPC , LC.zFPC]; % Focal Plane Center
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
        Xcenters = linspace(LC.FieldSize*(1/LC.res-1),LC.FieldSize*(1-1/LC.res),LC.res)/2;
        Ycenters = linspace(LC.FieldSize*(1/LC.res-1),LC.FieldSize*(1-1/LC.res),LC.res)/2;
        fprintf('%.3g%% of input power ends up on the detector.\n',100*mean(mean(sum(MCoutput.Image,3)))*LC.FieldSize^2);
    else
        fprintf('%.3g%% of input power ends up on the detector.\n',100*sum(MCoutput.Image,3));
    end

    if LC.nTimeBins > 0
        timevector = (-1/2:(LC.nTimeBins+1/2))*(LC.tEnd-LC.tStart)/LC.nTimeBins + LC.tStart;
    end

    if LC.res > 1 && LC.nTimeBins > 0
        h_f = plotVolumetric(7,Xcenters,Ycenters,timevector,MCoutput.Image,'slicePositions',[1 1 0]);
        h_f.Name = 'Image';
        xlabel('X [cm]');
        ylabel('Y [cm]');
        zlabel('Time [s]');
        title({'Normalized time-resolved fluence rate in the image plane','at 1x magnification [W/cm^2/W.incident]'});
        fprintf('Time-resolved light collector data plotted. Note that first time bin includes all\n  photons at earlier times and last time bin includes all photons at later times.\n');
    elseif LC.res > 1
        if(~ishandle(7))
            h_f = figure(7);
            h_f.Position = [40 80 1100 650];
        else
            h_f = figure(7);
        end
        clf;
        h_f.Name = 'Image';
        imagesc(Xcenters,Ycenters,MCoutput.Image.');
        title('Normalized fluence rate in the image plane at 1x magnification [W/cm^2/W.incident]');
        axis xy;axis equal;axis tight;xlabel('X [cm]');ylabel('Y [cm]');
        set(gca,'FontSize',18);
        colormap(inferno);
        colorbar;
    elseif LC.nTimeBins > 0
        if(~ishandle(7))
            h_f = figure(7);
            h_f.Position = [40 80 1100 650];
        else
            h_f = figure(7);
        end
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
drawnow;
end

