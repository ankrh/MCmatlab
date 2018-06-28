function [temperatureSensor, timeVector] = simulateHeatDistribution(name)
%   
%   Simulates temperature evolution due to absorption of light and 
%   diffusion of heat through the tissue.
%
%   Define the power, the on- and off-duration of the light source,
%   the number of pulses as well as the starting temperature distribution.
%   Define the number of times the temperatures should be extracted from
%   the simulation and plotted/recorded. Each update takes about 50-500 ms
%   to render depending on the number of voxels.
%   A movie wil be generated and saved if makemovie is set to true.
%   Note that the 3D matrices are defined in xyz-coordinates (and not yxz).
%   
%
%   Input
%       name
%           the basename of the files as specified in makeTissue.m
%
%   Displays
%       Tissue cuboid
%       3D temperature evolution during illumination and diffusion
%       Temperature evolution for the individual temperature sensors (x-y-plot)
%       3D temperature distributions after illumination and after diffusion
%       Tissue cuboid showing any heat-induced tissue damage
%
%   Output
%       temperatureSensor
%           array containing the temperatures of the temperature sensors at
%           each timestep
%       timeVector
%           vector containing the times the temperatures were recorded
%
%       ./Data/[name]_heatSimOutput.mat
%           file containing 3D temperature distribution after illumination
%           and after diffusion as well as temperatureSensor, timeVector
%           and other input values
%
%       ./Data/[name]_heatSimOutput.mp4
%           movie file showing the temperature evolution. The tissue cuboid
%           is shown in the beginning of the video.
%
%   Requires
%       deleteDataFiles.m
%       calcdtmax.m
%       plotVolumetric.m
%       updateVolumetric.m
%       finiteElementHeatPropagator.mex (architecture specific)
%

%% Updates
%   2014-01: Rasmus L. Pedersen & Mathias Christensen, DTU Fotonik
%   2017-06: Anders K. Hansen & Dominik Marti, DTU Fotonik

%% Check for preexisting files
if(~deleteDataFiles(name)); return; end

%% Load data from makeTissue.m and MonteCarlo.m
load(['./Data/' name '.mat']);
load(['./Data/' name '_MCoutput.mat'],'F');

%% Define parameters (user-specified)
Winc             = 4; % [W] Incident pulse peak power (in case of infinite plane waves, only the power incident upon the cuboid's top surface)
onduration       = 0.001; % [s] Pulse on-duration
offduration      = 0.004; % [s] Pulse off-duration
Temp             = 37*ones(size(T)); % [deg C] Initial temperature distribution
n_pulses         = 5; % Number of consecutive pulses, each with an illumination phase and a diffusion phase. If simulating only illumination or only diffusion, use n_pulses = 1.

plotTempLimits   = [37 105]; % [deg C], the expected range of temperatures, used only for setting the color scale in the plot
extremeprecision = false; % false: normal time step resolution (dt = dtmax/2), true: high time step resolution (dt = dtmax/100)
n_updates_on     = 30; % Number of times data is extracted for plots during each illumination phase . Must be at least 1.
n_updates_off    = 120; % Number of times data is extracted for plots during each diffusion phase. Must be at least 1.
makemovie        = true; % flag to specify whether a movie of the temperature evolution should be output. Remember to set plotTempLimits appropriately.

%% Determine remaining parameters
R = 8.3144598; % Gas constant [J/mol/K]

if n_pulses ~= 1 && (onduration == 0 || offduration == 0)
    n_pulses = 1; % If either pulse on-duration or off-duration is 0, it doesn't make sense to talk of multiple pulses
    fprintf('\nWarning: Number of pulses changed to 1 since either on-duration or off-duration is 0.\n\n');
end

[nx,ny,nz] = size(T);
nT = length(tissueList); % Number of different tissues in simulation

dx = x(2) - x(1);
dy = y(2) - y(1);
dz = z(2) - z(1);

HC = dx*dy*dz*[tissueList.VHC]; % Heat capacity array. Element i is the HC of tissue i in the reduced tissue list.
TC = [tissueList.TC]; % Thermal conductivity array. Element i is the TC of tissue i in the reduced tissue list.

if(any([tissueList.A])) % If non-zero Arrhenius data exists, prepare to calculate tissue damage.
    calcDamage = true;
    A = [tissueList.A];
    E = [tissueList.E];
    Omega = zeros(size(T));
else
    calcDamage = false;
end

%% Calculate time step
if extremeprecision
    dtmax = calcdtmax(T,TC,HC,dx,dy,dz)/100;
else
    dtmax = calcdtmax(T,TC,HC,dx,dy,dz)/2;
end

if onduration ~= 0
    nt_on = max(n_updates_on,ceil(onduration/dtmax)); % Number of time steps with illumination
    dt = onduration/nt_on; % Time step size
    if offduration ~= 0
        nt_off = max(n_updates_off,ceil(offduration/dt)); % Number of time steps without illumination
    else
        nt_off = 0;
    end
else
    nt_on = 0;
    nt_off = max(n_updates_off,ceil(offduration/dtmax)); % Number of time steps without illumination
    dt = offduration/nt_off; % Time step size
end

nt_vec = unique(floor([linspace(0,nt_on,n_updates_on+1) linspace(nt_on + nt_off/n_updates_off,nt_on + nt_off,n_updates_off)])); % Vector of nt's on which the plots should be updated

timeVector = [0 , repmat(nt_vec(2:end),1,n_pulses)+repelem(0:n_pulses-1,length(nt_vec)-1)*nt_vec(end)]*dt;

%% Calculate proportionality between voxel-to-voxel temperature difference DeltaT and time step temperature change dT
TC_eff = 2*(TC'*TC)./(ones(length(TC))*diag(TC)+diag(TC)*ones(length(TC))); % Same as TC_eff(i,j) = 2*TC_red(i)*TC_red(j)/(TC_red(i)+TC_red(j)) but without for loops
TC_eff(isnan(TC_eff)) = 0; % Neighboring insulating voxels return NaN but should just be 0
dTperdeltaT = cat(3,dt/dx*dy*dz*TC_eff./(diag(HC)*ones(nT)),dt*dx/dy*dz*TC_eff./(diag(HC)*ones(nT)),dt*dx*dy/dz*TC_eff./(diag(HC)*ones(nT))); % Third dimension corresponds to the different directions of heat diffusion
clear TC_eff

%% Calculate temperature change due to absorbed heat per time step
mua_vec = [tissueList.mua];
dT_abs = mua_vec(T).*F*dt*dx*dy*dz*Winc./HC(T); % Temperature change from absorption per time step [deg C] (mua_vec(T).*F is a 3D matrix of normalized volumetric powers)
clear F

%% Prepare the temperature plot
if(~ishandle(9))
    heatsimFigure = figure(9);
    heatsimFigure.Position = [40 80 1100 650];
else
    heatsimFigure = figure(9);
end
clf;
heatsimFigure.Name = 'Temperature evolution';
plotVolumetric(x,y,z,Temp,'reverseZ');
h_title = title(['Temperature evolution, t = ' num2str(timeVector(1),'%#.2g') ' s']);
caxis(plotTempLimits); % User-defined color scale limits

%% Make plots to visualize tissue properties
if(~ishandle(1))
    tissueFigure = figure(1);
    tissueFigure.Position = [40 80 1100 650];
else
    tissueFigure = figure(1);
end
clf;
tissueFigure.Name = 'Tissue type illustration';
plotVolumetric(x,y,z,T,'reverseZTissueIllustration',tissueList);
title('Tissue type illustration');

datacursormode on;
dataCursorHandle = datacursormode(tissueFigure);

fprintf('Place one or more (shift+click) temperature sensor in the tissue,\n and hit (almost) any key to continue.\nYou can also resize the temperature figure now and position the slices,\n which is useful for movie generation.\n');
pause
cursorInfo = getCursorInfo(dataCursorHandle);
numTemperatureSensors = length(cursorInfo);

if numTemperatureSensors
    for temperatureSensorIndex = numTemperatureSensors:-1:1
        dataCursorPosition = cursorInfo(temperatureSensorIndex).Position;
        [~, dCPz] = min(abs(z-dataCursorPosition(3)));
        [~, dCPy] = min(abs(y-dataCursorPosition(2)));
        [~, dCPx] = min(abs(x-dataCursorPosition(1)));
        temperatureSensorPosition(temperatureSensorIndex) = ...
            sub2ind(size(T),dCPx,dCPy,dCPz);
    end
    temperatureSensorPosition = fliplr(temperatureSensorPosition); % Reverse array so the order of the sensors is the same as the order they were defined in.
    temperatureSensor = NaN(numTemperatureSensors,length(timeVector));
    temperatureSensor(:,1) = Temp(temperatureSensorPosition);
    
    temperatureSensorTissues = {tissueList(T(temperatureSensorPosition)).name};
    
    if(~ishandle(11))
        temperatureSensorFigure = figure(11);
        temperatureSensorFigure.Position = [40 80 1100 650];
    else
        temperatureSensorFigure = figure(11);
    end
    clf;
    temperatureSensorFigure.Name = 'Temperature sensors';
    temperaturePlot = axes('XLim',timeVector([1, end]));
    xlabel('Time [sec]')
    ylabel('Temperature [deg C]')
    title('Temperature sensors')
    hold on
    temperatureSensorLinePlots = plot(temperaturePlot,timeVector,temperatureSensor);
    for i=1:numTemperatureSensors
        temperatureSensorLinePlots(i).YDataSource = sprintf('temperatureSensor(%d,:)',i);
    end
    legend(temperatureSensorTissues);
else
    temperatureSensor = [];
end

%% Simulate heat transfer
if numTemperatureSensors; figure(temperatureSensorFigure); end

fprintf('[nx,ny,nz]=[%d,%d,%d]. Number of pulses is %d.\nIllumination on for %d steps and off for %d steps in each pulse. Step size is %0.2e s.\n',nx,ny,nz,n_pulses,nt_on,nt_off,dt);

if(makemovie) % Make a temporary figure showing the tissue type illustration to put into the beginning of the movie
    tempFigure = figure;
    plotVolumetric(x,y,z,T,'reverseZTissueIllustration',tissueList);
    title('Tissue type illustration');
    tempFigure.Position = heatsimFigure.Position; % Size has to be the same as the temperature plot
    tempFigure.Children(6).Value = heatsimFigure.Children(6).Value; % Slices have to be set at the same positions
    tempFigure.Children(7).Value = heatsimFigure.Children(7).Value;
    tempFigure.Children(8).Value = heatsimFigure.Children(8).Value;
    callbackfunc = tempFigure.Children(6).Callback{1};
    vars = tempFigure.Children(6).Callback{2};
    feval(callbackfunc,tempFigure.Children(6),[],vars); % To update the display of the slices we have to call the callback function
    feval(callbackfunc,tempFigure.Children(7),[],vars);
    feval(callbackfunc,tempFigure.Children(8),[],vars);
    drawnow;
    movieframes(1) = getframe(tempFigure);
    delete(tempFigure);
    
    movieframes(2) = getframe(heatsimFigure);
end

tic
for j=1:n_pulses
    if nt_on
        fprintf(['Illuminating pulse #' num2str(j) '... \n' repmat('-',1,n_updates_on)]);
    else
        fprintf(['Diffusing heat... \n' repmat('-',1,n_updates_off)]);
    end
    drawnow;
    for i = 2:length(nt_vec)
        frameidx = i+(j-1)*(length(nt_vec)-1);
        Temp = finiteElementHeatPropagator(nt_vec(i)-nt_vec(i-1),Temp,T-1,dTperdeltaT,(nt_vec(i) <= nt_on)*dT_abs); % Arguments (nt,[[[Temp]]],[[[T]]],[[[dTperdeltaT]]],[[[dT_abs]]]). Contents of T have to be converted from Matlab's 1-based indexing to C's 0-based indexing.
        if calcDamage; Omega = Omega + (timeVector(frameidx)-timeVector(frameidx-1))*A(T).*exp(-E(T)./(R*(Temp + 273.15))); end;
        
        if numTemperatureSensors
            temperatureSensor(:,frameidx) = Temp(temperatureSensorPosition);
            refreshdata(temperatureSensorLinePlots,'caller');
        end
        
        fprintf(1,'\b');
        updateVolumetric(heatsimFigure,Temp);
        h_title.String = ['Temperature evolution, t = ' num2str(timeVector(frameidx),'%#.2g') ' s'];
        
        if(makemovie)
            movieframes(frameidx+1) = getframe(heatsimFigure);
        end
        
        if nt_vec(i) == nt_on
            if n_pulses == 1
                Temp_illum = Temp;
            end
            fprintf('\b\b Done\n');
            if nt_off
                fprintf(['Diffusing heat... \n' repmat('-',1,n_updates_off)]);
            end
        elseif nt_vec(i) == nt_on + nt_off
            fprintf('\b\b Done\n');
        end
    end
end
toc;

%% Plot and save results
figure(heatsimFigure)
save(['./Data/' name '_heatSimoutput.mat'],'temperatureSensor','timeVector','Winc','onduration','offduration','Temp','n_pulses');
if ~nt_on
    heatsimFigure.Name = 'Temperature after diffusion';
    title('Temperature after diffusion');
elseif ~nt_off
    heatsimFigure.Name = 'Temperature after illumination';
    title('Temperature after illumination');
elseif n_pulses == 1
    heatsimFigure.Name = 'Temperature after diffusion';
    title('Temperature after diffusion')
    if(~ishandle(10))
        illumFigure = figure(10);
        illumFigure.Position = [40 80 1100 650];
    else
        illumFigure = figure(10);
    end
    clf;
    illumFigure.Name = 'Temperature after illumination';
    plotVolumetric(x,y,z,Temp_illum,'reverseZ');
    title('Temperature after illumination');
    caxis(plotTempLimits); % User-defined color scale limits
    save(['./Data/' name '_heatSimoutput.mat'],'Temp_illum','-append');
else
    title(['Temperature after ' num2str(n_pulses) ' pulses'])
end

if calcDamage
    if(~ishandle(12))
        damageFigure = figure(12);
        damageFigure.Position = [40 80 1100 650];
    else
        damageFigure = figure(12);
    end
    damageFigure.Name = 'Tissue damage illustration';
    T_damage = T;
    T_damage(Omega > 1) = nT + 1;
    tissueList(nT + 1).name = 'damage';
    plotVolumetric(x,y,z,T_damage,'reverseZTissueIllustration',tissueList);
    title('Tissue damage illustration');
    save(['./Data/' name '_heatSimoutput.mat'],'Omega','-append');
end
fprintf('./Data/%s_heatSimoutput.mat saved\n',name);

if(makemovie)
    movieframes = [repmat(movieframes(1),1,30) movieframes(1:end) repmat(movieframes(end),1,30)];
    writerObj = VideoWriter(['./Data/' name '_heatSimoutput.mp4'],'MPEG-4');
    writerObj.Quality = 100;
    open(writerObj);
    warning('off','MATLAB:audiovideo:VideoWriter:mp4FramePadded');
    writeVideo(writerObj,movieframes);
    warning('on','MATLAB:audiovideo:VideoWriter:mp4FramePadded');
    close(writerObj);
    fprintf('./Data/%s_heatSimoutput.mp4 saved\n',name);
end

end
