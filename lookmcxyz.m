function lookmcxyz(name)
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

%% Updates
%   2014-08: Steven L. Jacques
%   2017-02: Steven L. Jacques
%   2017-06: Anders K. Hansen & Dominik Marti, DTU Fotonik
%   2018-04: Anders K. Hansen

load(['./Data/' name '.mat']);

%% Make tissue plot
if(~ishandle(1))
    h_f = figure(1);
    h_f.Position = [40 80 1100 650];
else
    h_f = figure(1);
end
h_f.Name = 'Tissue type illustration';
plotVolumetric(x,y,z,T,tissueList);
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
    load(['./Data/' name '_MCoutput.mat'],'F');
    dx = x(2)-x(1); dy = y(2)-y(1); dz = z(2)-z(1);
    
    %% Make fluence rate plot
    if(~ishandle(4))
        h_f = figure(4);
        h_f.Position = [40 80 1100 650];
    else
        h_f = figure(4);
    end
    h_f.Name = 'Normalized fluence rate';
    plotVolumetric(x,y,z,F);
    title('Normalized fluence rate (Intensity) [W/cm^2/W.incident] ')
    
    %% Make power absorption plot
    if(~ishandle(5))
        h_f = figure(5);
        h_f.Position = [40 80 1100 650];
    else
        h_f = figure(5);
    end
    h_f.Name = 'Normalized power absorption';
    mua_vec = [tissueList.mua];
    plotVolumetric(x,y,z,mua_vec(T).*F);
    title('Normalized absorbed power per unit volume [W/cm^3/W.incident] ')
    
    fprintf('\n%.2g%% of the input light was absorbed within the volume.\n',100*dx*dy*dz*sum(mua_vec(T(:)).*F(:).'));

    if(exist(['./Data/' name '_MCoutput_fluorescence.mat'],'file'))
        load(['./Data/' name '_MCoutput_fluorescence.mat'],'P','I_fluorescence');
        
        %% Remind the user what the input power was and plot emitter distribution
        fprintf('\nFluorescence was simulated for %.2g W of input excitation power.\n',P);
        
        fprintf('Out of this, %.2g W was absorbed within the volume.\n',dx*dy*dz*sum(mua_vec(T(:)).*P.*F(:).'));
        
        if(~ishandle(6))
            h_f = figure(6);
            h_f.Position = [40 80 1100 650];
        else
            h_f = figure(6);
        end
        h_f.Name = 'Fluorescence emitters';
        Y_vec = [tissueList.Y]; % The tissues' fluorescence power efficiencies
        sat_vec = [tissueList.sat]; % The tissues' fluorescence saturation fluence rates (intensity)
        FluorescenceEmitters = Y_vec(T).*mua_vec(T)*P.*F./(1 + P*F./sat_vec(T)); % [W/cm^3]
        plotVolumetric(x,y,z,FluorescenceEmitters);
        title('Fluorescence emitter distribution [W/cm^3] ')

        fprintf('Out of this, %.2g W was re-emitted as fluorescence.\n',dx*dy*dz*sum(FluorescenceEmitters(:)));
        
        %% Make fluence rate plot
        if(~ishandle(7))
            h_f = figure(7);
            h_f.Position = [40 80 1100 650];
        else
            h_f = figure(7);
        end
        h_f.Name = 'Fluorescence fluence rate';
        plotVolumetric(x,y,z,I_fluorescence);
        title('Fluorescence fluence rate (Intensity) [W/cm^2] ')
        
        %% Make power absorption plot
        if(~ishandle(8))
            h_f = figure(8);
            h_f.Position = [40 80 1100 650];
        else
            h_f = figure(8);
        end
        h_f.Name = 'Fluorescence power absorption';
        mua_vec = [tissueList_fluorescence.mua];
        plotVolumetric(x,y,z,mua_vec(T).*I_fluorescence);
        title('Absorbed fluorescence power per unit volume [W/cm^3] ')

        fprintf('Out of this, %.2g W was re-absorbed within the volume.\n\n',dx*dy*dz*sum(mua_vec(T(:)).*I_fluorescence(:).'));
    end
end

end
