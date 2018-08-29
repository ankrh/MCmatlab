function plotMCmatlabFluorescence(FMCinput,FMCoutput)
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

G = FMCinput.G;
mua_vec = [G.mediaProperties.mua];

%% Remind the user what the input power was and plot emitter distribution
P = FMCinput.Beam.P_excitation;
fprintf('\nFluorescence was simulated for %.2g W of input excitation power.\n',P);

fprintf('Out of this, %.3g W was absorbed within the cuboid.\n',G.dx*G.dy*G.dz*sum(sum(sum(mua_vec(G.M).*P.*FMCinput.MCoutput.F))));

Y_vec = [G.mediaProperties.Y]; % The media's fluorescence power efficiencies
sat_vec = [G.mediaProperties.sat]; % The media's fluorescence saturation fluence rates (intensity)
FluorescenceEmitters = Y_vec(G.M).*mua_vec(G.M)*P.*FMCinput.MCoutput.F./(1 + P*FMCinput.MCoutput.F./sat_vec(G.M)); % [W/cm^3]
h_f = plotVolumetric(8,G.x,G.y,G.z,FluorescenceEmitters,'MCmatlab_fromZero');
h_f.Name = 'Fluorescence emitters';
title('Fluorescence emitter distribution [W/cm^3] ')

fprintf('Out of this, %.3g W was re-emitted as fluorescence.\n',G.dx*G.dy*G.dz*sum(sum(sum(FluorescenceEmitters))));

%% Make power absorption plot
mua_vec = [G.mediaProperties_f.mua];
h_f = plotVolumetric(9,G.x,G.y,G.z,mua_vec(G.M).*FMCoutput.F,'MCmatlab_fromZero');
h_f.Name = 'Fluorescence power absorption';
title('Absorbed fluorescence power per unit volume [W/cm^3] ')

%% Make fluence rate plot
h_f = plotVolumetric(10,G.x,G.y,G.z,FMCoutput.F,'MCmatlab_fromZero');
h_f.Name = 'Fluorescence fluence rate';
title('Fluorescence fluence rate (Intensity) [W/cm^2] ')

fprintf('Out of this, %.3g W was re-absorbed within the cuboid.\n\n',G.dx*G.dy*G.dz*sum(sum(sum(mua_vec(G.M).*FMCoutput.F))));

if(FMCinput.useLightCollector)
    %% If there's a fluorescence light collector, show its orientation and the detected light
    h_f = plotVolumetric(11,G.x,G.y,G.z,G.M,'MCmatlab_GeometryIllustration',G.mediaProperties_f);
    h_f.Name = 'Fluorescence light collector illustration';
    title('Fluorescence light collector illustration');
    box on;grid on;grid minor;

    LC = FMCinput.LightCollector;
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
        fprintf('%.3g%% of fluorescence ends up on the detector.\n',100*mean(mean(FMCoutput.Image))*LC.FieldSize^2);
    else
        fprintf('%.3g%% of fluorescence ends up on the detector.\n',100*FMCoutput.Image);
    end

    if LC.res > 1
        if(~ishandle(12))
            h_f = figure(12);
            h_f.Position = [40 80 1100 650];
        else
            h_f = figure(12);
        end
        clf;
        h_f.Name = 'Fluorescence image';
        imagesc(Xcenters,Ycenters,FMCoutput.Image.');
        title('Fluence rate in the fluorescence image plane at 1x magnification [W/cm^2]');axis xy;axis equal;axis tight;xlabel('X [cm]');ylabel('Y [cm]');
        set(gca,'FontSize',18);
        colormap(inferno);
        colorbar;
    end
end
end
