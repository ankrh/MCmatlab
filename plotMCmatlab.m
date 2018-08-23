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
h_f = plotVolumetric(1,G.x,G.y,G.z,G.M,'MCmatlab_GeometryIllustration',G.mediaProperties);
h_f.Name = 'Geometry illustration';
title('Geometry illustration');

%% Make media properties plot
h_f = plotMediaProperties(2,G.mediaProperties);
h_f.Name = 'Media properties';

if(~isnan(G.wavelength_f))
    %% Make fluorescence media properties plot
    h_f = plotMediaProperties(3,G.mediaProperties_f);
    h_f.Name = 'Fluorescence media properties';
end

if(exist(['./Data/' name '_MCoutput.mat'],'file'))
    load(['./Data/' name '_MCoutput.mat'],'MCoutput','MCinput');
    
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
        else
            LCC = FPC;
            detectoraperture = LC.diam_LC/2*(cos(linspace(0,2*pi,100).')*Xvec + sin(linspace(0,2*pi,100).')*Yvec) + LCC;
            h2 = line(detectoraperture(:,1),detectoraperture(:,2),detectoraperture(:,3),'Color','r','LineWidth',2);
            legend(h2,'Fiber aperture','Location','northeast');
        end
        
        if LC.res_LC > 1
            Xcenters = linspace(LC.FieldSize_LC*(1/LC.res_LC-1),LC.FieldSize_LC*(1-1/LC.res_LC),LC.res_LC)/2;
            Ycenters = linspace(LC.FieldSize_LC*(1/LC.res_LC-1),LC.FieldSize_LC*(1-1/LC.res_LC),LC.res_LC)/2;
            fprintf('%.3g%% of input power ends up on the detector.\n',100*mean(mean(sum(MCoutput.Image,3)))*LC.FieldSize_LC^2);
        else
            fprintf('%.3g%% of input power ends up on the detector.\n',100*sum(MCoutput.Image,3));
        end
        
        if LC.nTimeBins_LC > 0
            timevector = (-1/2:(LC.nTimeBins_LC+1/2))*(LC.tEnd_LC-LC.tStart_LC)/LC.nTimeBins_LC + LC.tStart_LC;
        end
        
        if LC.res_LC > 1 && LC.nTimeBins_LC > 0
            h_f = plotVolumetric(7,Xcenters,Ycenters,timevector,MCoutput.Image,'slicePositions',[1 1 0]);
            h_f.Name = 'Image';
            xlabel('X [cm]');
            ylabel('Y [cm]');
            zlabel('Time [s]');
            title({'Normalized time-resolved fluence rate in the image plane','at 1x magnification [W/cm^2/W.incident]'});
            fprintf('Time-resolved light collector data plotted. Note that first time bin includes all\n  photons at earlier times and last time bin includes all photons at later times.\n');
        elseif LC.res_LC > 1
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
            colormap(GPBGYRcolormap);
            colorbar;
        elseif LC.nTimeBins_LC > 0
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

    if(exist(['./Data/' name '_MCoutput_fluorescence.mat'],'file'))
        load(['./Data/' name '_MCoutput_fluorescence.mat'],'MCoutput_f','MCinput_f');
        
        %% Remind the user what the input power was and plot emitter distribution
        P = MCinput_f.Beam.P_excitation;
        fprintf('\nFluorescence was simulated for %.2g W of input excitation power.\n',P);
        
        fprintf('Out of this, %.3g W was absorbed within the cuboid.\n',G.dx*G.dy*G.dz*sum(sum(sum(mua_vec(G.M).*P.*MCoutput.F))));
        
        Y_vec = [G.mediaProperties.Y]; % The media's fluorescence power efficiencies
        sat_vec = [G.mediaProperties.sat]; % The media's fluorescence saturation fluence rates (intensity)
        FluorescenceEmitters = Y_vec(G.M).*mua_vec(G.M)*P.*MCoutput.F./(1 + P*MCoutput.F./sat_vec(G.M)); % [W/cm^3]
        h_f = plotVolumetric(8,G.x,G.y,G.z,FluorescenceEmitters,'MCmatlab_fromZero');
        h_f.Name = 'Fluorescence emitters';
        title('Fluorescence emitter distribution [W/cm^3] ')

        fprintf('Out of this, %.3g W was re-emitted as fluorescence.\n',G.dx*G.dy*G.dz*sum(sum(sum(FluorescenceEmitters))));
        
        %% Make power absorption plot
        mua_vec = [G.mediaProperties_f.mua];
        h_f = plotVolumetric(9,G.x,G.y,G.z,mua_vec(G.M).*MCoutput_f.F,'MCmatlab_fromZero');
        h_f.Name = 'Fluorescence power absorption';
        title('Absorbed fluorescence power per unit volume [W/cm^3] ')

        %% Make fluence rate plot
        h_f = plotVolumetric(10,G.x,G.y,G.z,MCoutput_f.F,'MCmatlab_fromZero');
        h_f.Name = 'Fluorescence fluence rate';
        title('Fluorescence fluence rate (Intensity) [W/cm^2] ')
        
        fprintf('Out of this, %.3g W was re-absorbed within the cuboid.\n\n',G.dx*G.dy*G.dz*sum(sum(sum(mua_vec(G.M).*MCoutput_f.F))));
        
        if(MCinput_f.useLightCollector)
            %% If there's a fluorescence light collector, show its orientation and the detected light
            h_f = plotVolumetric(11,G.x,G.y,G.z,G.M,'MCmatlab_GeometryIllustration',G.mediaProperties_f);
            h_f.Name = 'Fluorescence light collector illustration';
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
            else
                LCC = FPC;
                detectoraperture = LC_f.diam_LC/2*(cos(linspace(0,2*pi,100).')*Xvec + sin(linspace(0,2*pi,100).')*Yvec) + LCC;
                h2 = line(detectoraperture(:,1),detectoraperture(:,2),detectoraperture(:,3),'Color','r','LineWidth',2);
                legend(h2,'Fiber aperture','Location','northeast');
            end
            
            if LC_f.res_LC > 1
                Xcenters_f = linspace(LC_f.FieldSize_LC*(1/LC_f.res_LC-1),LC_f.FieldSize_LC*(1-1/LC_f.res_LC),LC_f.res_LC)/2;
                Ycenters_f = linspace(LC_f.FieldSize_LC*(1/LC_f.res_LC-1),LC_f.FieldSize_LC*(1-1/LC_f.res_LC),LC_f.res_LC)/2;
                fprintf('%.3g%% of fluorescence ends up on the detector.\n',100*mean(mean(MCoutput_f.Image))*LC_f.FieldSize_LC^2);
            else
                fprintf('%.3g%% of fluorescence ends up on the detector.\n',100*MCoutput_f.Image);
            end
            
            if LC_f.res_LC > 1
                if(~ishandle(12))
                    h_f = figure(12);
                    h_f.Position = [40 80 1100 650];
                else
                    h_f = figure(12);
                end
                clf;
                h_f.Name = 'Fluorescence image';
                imagesc(Xcenters_f,Ycenters_f,MCoutput_f.Image.');
                title('Fluence rate in the fluorescence image plane at 1x magnification [W/cm^2]');axis xy;axis equal;axis tight;xlabel('X [cm]');ylabel('Y [cm]');
                set(gca,'FontSize',18);
                colormap(GPBGYRcolormap);
                colorbar;
            end
        end
    end
end
end
