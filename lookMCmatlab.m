function lookMCmatlab(name)
%
%   Displays the tissue cuboid, an overview over the known optical, thermal
%   and flourescence properties of the tissue types and the output of any
%   completed Monte Carlo simulations (excitation and/or fluorescence).
%
%   Input
%       name
%           the basename of the files as specified in makeTissue
%
%   Displays
%       Tissue cuboid
%       Tissue optical, thermal and fluorescence properties
%   If Monte Carlo output data exists, displays
%       Fluence rate
%       Absorbed power
%   And, if fluorescence Monte Carlo output data exists, displays
%       Tissue optical and thermal properties at the fluorescence wavelength
%       Distribution of fluorescence emitters
%       Fluorescence fluence rate
%       Absorbed fluorescence power
%
%   Requires
%       plotVolumetric.m
%       plotTissueProperties.m
%

%% Acknowledgement
%   This function was inspired by lookmcxyz of the mcxyz MC program hosted at omlc.org

%% Load tissue definition
load(['./Data/' name '.mat']);

%% Make tissue plot
if(~ishandle(1))
    h_f = figure(1);
    h_f.Position = [40 80 1100 650];
else
    h_f = figure(1);
end
h_f.Name = 'Tissue type illustration';
plotVolumetric(x,y,z,T,'MCmatlab_TissueIllustration',tissueList);
title('Tissue type illustration');

%% Make tissue properties plot
if(~ishandle(2))
    h_f = figure(2);
    h_f.Position = [40 80 1100 650];
else
    h_f = figure(2);
end
h_f.Name = 'Tissue properties';
plotTissueProperties(tissueList);

if(exist('tissueList_fluorescence','var'))
    %% Make fluorescence tissue properties plot
    if(~ishandle(3))
        h_f = figure(3);
        h_f.Position = [40 80 1100 650];
    else
        h_f = figure(3);
    end
    h_f.Name = 'Fluorescence tissue properties';
    plotTissueProperties(tissueList_fluorescence);
end

if(exist(['./Data/' name '_MCoutput.mat'],'file'))
    load(['./Data/' name '_MCoutput.mat'],'MCoutput','MCinput');
    dx = x(2)-x(1); dy = y(2)-y(1); dz = z(2)-z(1);
    
    %% Make fluence rate plot
    if(~ishandle(4))
        h_f = figure(4);
        h_f.Position = [40 80 1100 650];
    else
        h_f = figure(4);
    end
    h_f.Name = 'Normalized fluence rate';
    plotVolumetric(x,y,z,MCoutput.F,'MCmatlab_fromZero');
    title('Normalized fluence rate (Intensity) [W/cm^2/W.incident] ')
    
    %% Make power absorption plot
    if(~ishandle(8))
        h_f = figure(8);
        h_f.Position = [40 80 1100 650];
    else
        h_f = figure(8);
    end
    h_f.Name = 'Normalized power absorption';
    mua_vec = [tissueList.mua];
    plotVolumetric(x,y,z,mua_vec(T).*MCoutput.F,'MCmatlab_fromZero');
    title('Normalized absorbed power per unit volume [W/cm^3/W.incident] ')
    
    fprintf('\n%.3g%% of the input light was absorbed within the volume.\n',100*dx*dy*dz*sum(sum(sum(mua_vec(T).*MCoutput.F))));

    if(MCinput.useLightCollector)
        %% If there's a light collector, show its orientation and the detected light
        if(~ishandle(7))
            h_f = figure(7);
            h_f.Position = [40 80 1100 650];
        else
            h_f = figure(7);
        end
        h_f.Name = 'Light collector illustration';
        plotVolumetric(x,y,z,T,'MCmatlab_TissueIllustration',tissueList);
        title('Light collector illustration');
        box on;grid on;grid minor;

        arrowlength = sqrt((nx*dx)^2+(ny*dy)^2+(nz*dz)^2)/5;
        Zvec = [sin(MCinput.theta_LC)*cos(MCinput.phi_LC) , sin(MCinput.theta_LC)*sin(MCinput.phi_LC) , cos(MCinput.theta_LC)];
        Xvec = [sin(MCinput.phi_LC) , -cos(MCinput.phi_LC) , 0];
        Yvec = cross(Zvec,Xvec);
        FPC = [MCinput.xFPC_LC , MCinput.yFPC_LC , MCinput.zFPC_LC]; % Focal Plane Center
        FPC_X = FPC + arrowlength*Xvec;
        line([FPC(1) FPC_X(1)],[FPC(2) FPC_X(2)],[FPC(3) FPC_X(3)],'Linewidth',2,'Color','r')
        text(FPC_X(1),FPC_X(2),FPC_X(3),'X','HorizontalAlignment','center','FontSize',18)
        FPC_Y = FPC + arrowlength*Yvec;
        line([FPC(1) FPC_Y(1)],[FPC(2) FPC_Y(2)],[FPC(3) FPC_Y(3)],'Linewidth',2,'Color','r')
        text(FPC_Y(1),FPC_Y(2),FPC_Y(3),'Y','HorizontalAlignment','center','FontSize',18)
        
        if isfinite(MCinput.f_LC)
            fieldperimeter = MCinput.FSorNA_LC/2*(cos(linspace(0,2*pi,100).')*Xvec + sin(linspace(0,2*pi,100).')*Yvec) + FPC;
            h1 = line(fieldperimeter(:,1),fieldperimeter(:,2),fieldperimeter(:,3),'Color','b','LineWidth',2);
            LCC = FPC - Zvec*MCinput.f_LC; % Light Collector Center
            detectoraperture = MCinput.diam_LC/2*(cos(linspace(0,2*pi,100).')*Xvec + sin(linspace(0,2*pi,100).')*Yvec) + LCC;
            h2 = line(detectoraperture(:,1),detectoraperture(:,2),detectoraperture(:,3),'Color','r','LineWidth',2);
            legend([h1 h2],'Imaged area','Lens aperture','Location','northeast');
        
            if(~ishandle(5))
                h_f = figure(5);
                h_f.Position = [40 80 1100 650];
            else
                h_f = figure(5);
            end
            clf;
            h_f.Name = 'Image plane';
            imagesc([-MCinput.FSorNA_LC MCinput.FSorNA_LC]/2,[-MCinput.FSorNA_LC MCinput.FSorNA_LC]/2,MCoutput.ImP.');
            title('Fluence rate in the image plane at 1x magnification [W/cm^2/W.incident]');axis xy;axis equal;axis tight;xlabel('X [cm]');ylabel('Y [cm]');
            set(gca,'FontSize',18);
            colormap(GPBGYRcolormap);
            colorbar;
        else
            LCC = FPC;
            detectoraperture = MCinput.diam_LC/2*(cos(linspace(0,2*pi,100).')*Xvec + sin(linspace(0,2*pi,100).')*Yvec) + LCC;
            h2 = line(detectoraperture(:,1),detectoraperture(:,2),detectoraperture(:,3),'Color','r','LineWidth',2);
            legend(h2,'Fiber aperture','Location','northeast');
        end
        
        if(~ishandle(6))
            h_f = figure(6);
            h_f.Position = [40 80 1100 650];
        else
            h_f = figure(6);
        end
        clf;
        h_f.Name = 'Light collector plane';
        imagesc([-MCinput.diam_LC MCinput.diam_LC]/2,[-MCinput.diam_LC MCinput.diam_LC]/2,MCoutput.LCP.');
        title('Light in the light collector plane that ends up detected [W/cm^2/W.incident]');axis xy;axis equal;axis tight;xlabel('X [cm]');ylabel('Y [cm]');
        set(gca,'FontSize',18);
        colormap(GPBGYRcolormap);
        colorbar;

        fprintf('%.3g%% of input power ends up on the detector.\n',100*mean(mean(MCoutput.LCP))*MCinput.diam_LC^2);
    end

    if(exist(['./Data/' name '_MCoutput_fluorescence.mat'],'file'))
        load(['./Data/' name '_MCoutput_fluorescence.mat'],'P','MCoutput_fluorescence','MCinput_fluorescence');
        
        %% Remind the user what the input power was and plot emitter distribution
        fprintf('\nFluorescence was simulated for %.2g W of input excitation power.\n',P);
        
        fprintf('Out of this, %.3g W was absorbed within the volume.\n',dx*dy*dz*sum(sum(sum(mua_vec(T).*P.*MCoutput.F))));
        
        if(~ishandle(9))
            h_f = figure(9);
            h_f.Position = [40 80 1100 650];
        else
            h_f = figure(9);
        end
        h_f.Name = 'Fluorescence emitters';
        Y_vec = [tissueList.Y]; % The tissues' fluorescence power efficiencies
        sat_vec = [tissueList.sat]; % The tissues' fluorescence saturation fluence rates (intensity)
        FluorescenceEmitters = Y_vec(T).*mua_vec(T)*P.*MCoutput.F./(1 + P*MCoutput.F./sat_vec(T)); % [W/cm^3]
        plotVolumetric(x,y,z,FluorescenceEmitters,'MCmatlab_fromZero');
        title('Fluorescence emitter distribution [W/cm^3] ')

        fprintf('Out of this, %.3g W was re-emitted as fluorescence.\n',dx*dy*dz*sum(sum(sum(FluorescenceEmitters))));
        
        %% Make fluence rate plot
        if(~ishandle(10))
            h_f = figure(10);
            h_f.Position = [40 80 1100 650];
        else
            h_f = figure(10);
        end
        h_f.Name = 'Fluorescence fluence rate';
        plotVolumetric(x,y,z,MCoutput_fluorescence.F,'MCmatlab_fromZero');
        title('Fluorescence fluence rate (Intensity) [W/cm^2] ')
        
        %% Make power absorption plot
        if(~ishandle(11))
            h_f = figure(11);
            h_f.Position = [40 80 1100 650];
        else
            h_f = figure(11);
        end
        h_f.Name = 'Fluorescence power absorption';
        mua_vec = [tissueList_fluorescence.mua];
        plotVolumetric(x,y,z,mua_vec(T).*MCoutput_fluorescence.F,'MCmatlab_fromZero');
        title('Absorbed fluorescence power per unit volume [W/cm^3] ')

        fprintf('Out of this, %.3g W was re-absorbed within the volume.\n\n',dx*dy*dz*sum(sum(sum(mua_vec(T).*MCoutput_fluorescence.F))));
        
        if(MCinput_fluorescence.useLightCollector)
            %% If there's a fluorescence light collector, show its orientation and the detected light
            if(~ishandle(12))
                h_f = figure(12);
                h_f.Position = [40 80 1100 650];
            else
                h_f = figure(12);
            end
            h_f.Name = 'Fluorescence light collector illustration';
            plotVolumetric(x,y,z,T,'MCmatlab_TissueIllustration',tissueList_fluorescence);
            title('Fluorescence light collector illustration');
            box on;grid on;grid minor;

            arrowlength = sqrt((nx*dx)^2+(ny*dy)^2+(nz*dz)^2)/5;
            Zvec = [sin(MCinput_fluorescence.theta_LC)*cos(MCinput_fluorescence.phi_LC) , sin(MCinput_fluorescence.theta_LC)*sin(MCinput_fluorescence.phi_LC) , cos(MCinput_fluorescence.theta_LC)];
            Xvec = [sin(MCinput_fluorescence.phi_LC) , -cos(MCinput_fluorescence.phi_LC) , 0];
            Yvec = cross(Zvec,Xvec);
            FPC = [MCinput_fluorescence.xFPC_LC , MCinput_fluorescence.yFPC_LC , MCinput_fluorescence.zFPC_LC]; % Focal Plane Center
            FPC_X = FPC + arrowlength*Xvec;
            line([FPC(1) FPC_X(1)],[FPC(2) FPC_X(2)],[FPC(3) FPC_X(3)],'Linewidth',2,'Color','r')
            text(FPC_X(1),FPC_X(2),FPC_X(3),'X','HorizontalAlignment','center','FontSize',18)
            FPC_Y = FPC + arrowlength*Yvec;
            line([FPC(1) FPC_Y(1)],[FPC(2) FPC_Y(2)],[FPC(3) FPC_Y(3)],'Linewidth',2,'Color','r')
            text(FPC_Y(1),FPC_Y(2),FPC_Y(3),'Y','HorizontalAlignment','center','FontSize',18)

            if isfinite(MCinput_fluorescence.f_LC)
                fieldperimeter = MCinput_fluorescence.FSorNA_LC/2*(cos(linspace(0,2*pi,100).')*Xvec + sin(linspace(0,2*pi,100).')*Yvec) + FPC;
                h1 = line(fieldperimeter(:,1),fieldperimeter(:,2),fieldperimeter(:,3),'Color','b','LineWidth',2);
                LCC = FPC - Zvec*MCinput_fluorescence.f_LC; % Light Collector Center
                detectoraperture = MCinput_fluorescence.diam_LC/2*(cos(linspace(0,2*pi,100).')*Xvec + sin(linspace(0,2*pi,100).')*Yvec) + LCC;
                h2 = line(detectoraperture(:,1),detectoraperture(:,2),detectoraperture(:,3),'Color','r','LineWidth',2);
                legend([h1 h2],'Imaged area','Lens aperture','Location','northeast');

                if(~ishandle(13))
                    h_f = figure(13);
                    h_f.Position = [40 80 1100 650];
                else
                    h_f = figure(13);
                end
                clf;
                h_f.Name = 'Fluorescence image plane';
                imagesc([-MCinput_fluorescence.FSorNA_LC MCinput_fluorescence.FSorNA_LC]/2,[-MCinput_fluorescence.FSorNA_LC MCinput_fluorescence.FSorNA_LC]/2,MCoutput_fluorescence.ImP.');
                title('Detected fluorescence light in the image plane [W/cm^2]');axis xy;axis equal;axis tight;xlabel('X [cm]');ylabel('Y [cm]');
                set(gca,'FontSize',18);
                colormap(GPBGYRcolormap);
                colorbar;
            else
                LCC = FPC;
                detectoraperture = MCinput_fluorescence.diam_LC/2*(cos(linspace(0,2*pi,100).')*Xvec + sin(linspace(0,2*pi,100).')*Yvec) + LCC;
                h2 = line(detectoraperture(:,1),detectoraperture(:,2),detectoraperture(:,3),'Color','r','LineWidth',2);
                legend(h2,'Fiber aperture','Location','northeast');
            end

            if(~ishandle(14))
                h_f = figure(14);
                h_f.Position = [40 80 1100 650];
            else
                h_f = figure(14);
            end
            clf;
            h_f.Name = 'Fluorescence light collector plane';
            imagesc([-MCinput_fluorescence.diam_LC MCinput_fluorescence.diam_LC]/2,[-MCinput_fluorescence.diam_LC MCinput_fluorescence.diam_LC]/2,MCoutput_fluorescence.LCP.');
            title('Detected fluorescence light in the light collector plane [W/cm^2]');axis xy;axis equal;axis tight;xlabel('X [cm]');ylabel('Y [cm]');
            set(gca,'FontSize',18);
            colormap(GPBGYRcolormap);
            colorbar;
            
            fprintf('%.3g W of fluorescence ends up on the detector.\n',mean(mean(MCoutput_fluorescence.LCP))*MCinput_fluorescence.diam_LC^2);
        end
    end
end
end
