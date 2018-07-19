function plotMCmatlab(name)
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

%% Load geometry definition and media properties
load(['./Data/' name '.mat']);

%% Make geometry plot
if(~ishandle(1))
    h_f = figure(1);
    h_f.Position = [40 80 1100 650];
else
    h_f = figure(1);
end
h_f.Name = 'Geometry illustration';
plotVolumetric(G.x,G.y,G.z,G.M,'MCmatlab_GeometryIllustration',G.mediaProperties);
title('Geometry illustration');

%% Make media properties plot
if(~ishandle(2))
    h_f = figure(2);
    h_f.Position = [40 80 1100 650];
else
    h_f = figure(2);
end
h_f.Name = 'Media properties';
plotMediaProperties(G.mediaProperties);

if(~isnan(G.wavelength_f))
    %% Make fluorescence media properties plot
    if(~ishandle(3))
        h_f = figure(3);
        h_f.Position = [40 80 1100 650];
    else
        h_f = figure(3);
    end
    h_f.Name = 'Fluorescence media properties';
    plotMediaProperties(G.mediaProperties_f);
end

if(exist(['./Data/' name '_MCoutput.mat'],'file'))
    load(['./Data/' name '_MCoutput.mat'],'MCoutput','MCinput');
    
    %% Make power absorption plot
    if(~ishandle(4))
        h_f = figure(4);
        h_f.Position = [40 80 1100 650];
    else
        h_f = figure(4);
    end
    h_f.Name = 'Normalized power absorption';
    mua_vec = [G.mediaProperties.mua];
    plotVolumetric(G.x,G.y,G.z,mua_vec(G.M).*MCoutput.F,'MCmatlab_fromZero');
    title('Normalized absorbed power per unit volume [W/cm^3/W.incident] ')
    
    %% Make fluence rate plot
    if(~ishandle(5))
        h_f = figure(5);
        h_f.Position = [40 80 1100 650];
    else
        h_f = figure(5);
    end
    h_f.Name = 'Normalized fluence rate';
    plotVolumetric(G.x,G.y,G.z,MCoutput.F,'MCmatlab_fromZero');
    title('Normalized fluence rate (Intensity) [W/cm^2/W.incident] ')
    
    fprintf('\n%.3g%% of the input light was absorbed within the cuboid.\n',100*G.dx*G.dy*G.dz*sum(sum(sum(mua_vec(G.M).*MCoutput.F))));

    if(MCinput.useLightCollector)
        %% If there's a light collector, show its orientation and the detected light
        if(~ishandle(6))
            h_f = figure(6);
            h_f.Position = [40 80 1100 650];
        else
            h_f = figure(6);
        end
        h_f.Name = 'Light collector illustration';
        plotVolumetric(G.x,G.y,G.z,G.M,'MCmatlab_GeometryIllustration',G.mediaProperties);
        title('Light collector illustration');
        box on;grid on;grid minor;

        LC = MCinput.LightCollector;
        arrowlength = sqrt((G.nx*G.dx)^2+(G.ny*G.dy)^2+(G.nz*G.dz)^2)/5;
        Zvec = [sin(LC.theta_LC)*cos(LC.phi_LC) , sin(LC.theta_LC)*sin(LC.phi_LC) , cos(LC.theta_LC)];
        Xvec = [sin(LC.phi_LC) , -cos(LC.phi_LC) , 0];
        Yvec = cross(Zvec,Xvec);
        FPC = [LC.xFPC_LC , LC.yFPC_LC , LC.zFPC_LC]; % Focal Plane Center
        FPC_X = FPC + arrowlength*Xvec;
        line([FPC(1) FPC_X(1)],[FPC(2) FPC_X(2)],[FPC(3) FPC_X(3)],'Linewidth',2,'Color','r')
        text(FPC_X(1),FPC_X(2),FPC_X(3),'X','HorizontalAlignment','center','FontSize',18)
        FPC_Y = FPC + arrowlength*Yvec;
        line([FPC(1) FPC_Y(1)],[FPC(2) FPC_Y(2)],[FPC(3) FPC_Y(3)],'Linewidth',2,'Color','r')
        text(FPC_Y(1),FPC_Y(2),FPC_Y(3),'Y','HorizontalAlignment','center','FontSize',18)
        
        if isfinite(LC.f_LC)
            fieldperimeter = LC.FieldSize_LC/2*(cos(linspace(0,2*pi,100).')*Xvec + sin(linspace(0,2*pi,100).')*Yvec) + FPC;
            h1 = line(fieldperimeter(:,1),fieldperimeter(:,2),fieldperimeter(:,3),'Color','b','LineWidth',2);
            LCC = FPC - Zvec*LC.f_LC; % Light Collector Center
            detectoraperture = LC.diam_LC/2*(cos(linspace(0,2*pi,100).')*Xvec + sin(linspace(0,2*pi,100).')*Yvec) + LCC;
            h2 = line(detectoraperture(:,1),detectoraperture(:,2),detectoraperture(:,3),'Color','r','LineWidth',2);
            legend([h1 h2],'Imaged area','Lens aperture','Location','northeast');

            if LC.resX_LC*LC.resY_LC > 1
                if(~ishandle(7))
                    h_f = figure(7);
                    h_f.Position = [40 80 1100 650];
                else
                    h_f = figure(7);
                end
                clf;
                h_f.Name = 'Image';
                imagesc([-LC.FieldSize_LC LC.FieldSize_LC]/2,[-LC.FieldSize_LC LC.FieldSize_LC]/2,MCoutput.Image.');
                title('Normalized fluence rate in the image plane at 1x magnification [W/cm^2/W.incident]');axis xy;axis equal;axis tight;xlabel('X [cm]');ylabel('Y [cm]');
                set(gca,'FontSize',18);
                colormap(GPBGYRcolormap);
                colorbar;
            end
        else
            LCC = FPC;
            detectoraperture = LC.diam_LC/2*(cos(linspace(0,2*pi,100).')*Xvec + sin(linspace(0,2*pi,100).')*Yvec) + LCC;
            h2 = line(detectoraperture(:,1),detectoraperture(:,2),detectoraperture(:,3),'Color','r','LineWidth',2);
            legend(h2,'Fiber aperture','Location','northeast');
        end

        if length(MCoutput.Image) > 1
            fprintf('%.3g%% of input power ends up on the detector.\n',100*mean(mean(MCoutput.Image))*LC.FieldSize_LC^2);
        else
            fprintf('%.3g%% of input power ends up on the detector.\n',100*MCoutput.Image);
        end
    end

    if(exist(['./Data/' name '_MCoutput_fluorescence.mat'],'file'))
        load(['./Data/' name '_MCoutput_fluorescence.mat'],'MCoutput_f','MCinput_f');
        
        %% Remind the user what the input power was and plot emitter distribution
        P = MCinput_f.Beam.P_excitation;
        fprintf('\nFluorescence was simulated for %.2g W of input excitation power.\n',P);
        
        fprintf('Out of this, %.3g W was absorbed within the cuboid.\n',G.dx*G.dy*G.dz*sum(sum(sum(mua_vec(G.M).*P.*MCoutput.F))));
        
        if(~ishandle(8))
            h_f = figure(8);
            h_f.Position = [40 80 1100 650];
        else
            h_f = figure(8);
        end
        h_f.Name = 'Fluorescence emitters';
        Y_vec = [G.mediaProperties.Y]; % The media's fluorescence power efficiencies
        sat_vec = [G.mediaProperties.sat]; % The media's fluorescence saturation fluence rates (intensity)
        FluorescenceEmitters = Y_vec(G.M).*mua_vec(G.M)*P.*MCoutput.F./(1 + P*MCoutput.F./sat_vec(G.M)); % [W/cm^3]
        plotVolumetric(G.x,G.y,G.z,FluorescenceEmitters,'MCmatlab_fromZero');
        title('Fluorescence emitter distribution [W/cm^3] ')

        fprintf('Out of this, %.3g W was re-emitted as fluorescence.\n',G.dx*G.dy*G.dz*sum(sum(sum(FluorescenceEmitters))));
        
        %% Make power absorption plot
        if(~ishandle(9))
            h_f = figure(9);
            h_f.Position = [40 80 1100 650];
        else
            h_f = figure(9);
        end
        h_f.Name = 'Fluorescence power absorption';
        mua_vec = [G.mediaProperties_f.mua];
        plotVolumetric(G.x,G.y,G.z,mua_vec(G.M).*MCoutput_f.F,'MCmatlab_fromZero');
        title('Absorbed fluorescence power per unit volume [W/cm^3] ')

        %% Make fluence rate plot
        if(~ishandle(10))
            h_f = figure(10);
            h_f.Position = [40 80 1100 650];
        else
            h_f = figure(10);
        end
        h_f.Name = 'Fluorescence fluence rate';
        plotVolumetric(G.x,G.y,G.z,MCoutput_f.F,'MCmatlab_fromZero');
        title('Fluorescence fluence rate (Intensity) [W/cm^2] ')
        
        fprintf('Out of this, %.3g W was re-absorbed within the cuboid.\n\n',G.dx*G.dy*G.dz*sum(sum(sum(mua_vec(G.M).*MCoutput_f.F))));
        
        if(MCinput_f.useLightCollector)
            %% If there's a fluorescence light collector, show its orientation and the detected light
            if(~ishandle(11))
                h_f = figure(11);
                h_f.Position = [40 80 1100 650];
            else
                h_f = figure(11);
            end
            h_f.Name = 'Fluorescence light collector illustration';
            plotVolumetric(G.x,G.y,G.z,G.M,'MCmatlab_GeometryIllustration',G.mediaProperties_f);
            title('Fluorescence light collector illustration');
            box on;grid on;grid minor;

            LC_f = MCinput_f.LightCollector;
            arrowlength = sqrt((G.nx*G.dx)^2+(G.ny*G.dy)^2+(G.nz*G.dz)^2)/5;
            Zvec = [sin(LC_f.theta_LC)*cos(LC_f.phi_LC) , sin(LC_f.theta_LC)*sin(LC_f.phi_LC) , cos(LC_f.theta_LC)];
            Xvec = [sin(LC_f.phi_LC) , -cos(LC_f.phi_LC) , 0];
            Yvec = cross(Zvec,Xvec);
            FPC = [LC_f.xFPC_LC , LC_f.yFPC_LC , LC_f.zFPC_LC]; % Focal Plane Center
            FPC_X = FPC + arrowlength*Xvec;
            line([FPC(1) FPC_X(1)],[FPC(2) FPC_X(2)],[FPC(3) FPC_X(3)],'Linewidth',2,'Color','r')
            text(FPC_X(1),FPC_X(2),FPC_X(3),'X','HorizontalAlignment','center','FontSize',18)
            FPC_Y = FPC + arrowlength*Yvec;
            line([FPC(1) FPC_Y(1)],[FPC(2) FPC_Y(2)],[FPC(3) FPC_Y(3)],'Linewidth',2,'Color','r')
            text(FPC_Y(1),FPC_Y(2),FPC_Y(3),'Y','HorizontalAlignment','center','FontSize',18)

            if isfinite(LC_f.f_LC)
                fieldperimeter = LC_f.FieldSize_LC/2*(cos(linspace(0,2*pi,100).')*Xvec + sin(linspace(0,2*pi,100).')*Yvec) + FPC;
                h1 = line(fieldperimeter(:,1),fieldperimeter(:,2),fieldperimeter(:,3),'Color','b','LineWidth',2);
                LCC = FPC - Zvec*LC_f.f_LC; % Light Collector Center
                detectoraperture = LC_f.diam_LC/2*(cos(linspace(0,2*pi,100).')*Xvec + sin(linspace(0,2*pi,100).')*Yvec) + LCC;
                h2 = line(detectoraperture(:,1),detectoraperture(:,2),detectoraperture(:,3),'Color','r','LineWidth',2);
                legend([h1 h2],'Imaged area','Lens aperture','Location','northeast');

                if LC_f.resX_LC*LC_f.resY_LC > 1
                    if(~ishandle(12))
                        h_f = figure(12);
                        h_f.Position = [40 80 1100 650];
                    else
                        h_f = figure(12);
                    end
                    clf;
                    h_f.Name = 'Fluorescence image';
                    imagesc([-LC_f.FieldSize_LC LC_f.FieldSize_LC]/2,[-LC_f.FieldSize_LC LC_f.FieldSize_LC]/2,MCoutput_f.Image.');
                    title('Fluence rate in the fluorescence image plane at 1x magnification [W/cm^2]');axis xy;axis equal;axis tight;xlabel('X [cm]');ylabel('Y [cm]');
                    set(gca,'FontSize',18);
                    colormap(GPBGYRcolormap);
                    colorbar;
                end
            else
                LCC = FPC;
                detectoraperture = LC_f.diam_LC/2*(cos(linspace(0,2*pi,100).')*Xvec + sin(linspace(0,2*pi,100).')*Yvec) + LCC;
                h2 = line(detectoraperture(:,1),detectoraperture(:,2),detectoraperture(:,3),'Color','r','LineWidth',2);
                legend(h2,'Fiber aperture','Location','northeast');
            end

            if length(MCoutput_f.Image) > 1
                fprintf('%.3g%% of fluorescence ends up on the detector.\n',100*mean(mean(MCoutput_f.Image))*LC_f.FieldSize_LC^2);
            else
                fprintf('%.3g%% of fluorescence ends up on the detector.\n',100*MCoutput_f.Image);
            end
        end
    end
end
end