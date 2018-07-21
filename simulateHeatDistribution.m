function [temperatureSensor, timeVector] = simulateHeatDistribution(name)
%   Created 2018 by Dominik Marti and Anders K. Hansen, DTU Fotonik
%   Work on this heat solver was started by Rasmus L. Pedersen & Mathias Christensen, DTU Fotonik
%   
%   Simulates temperature evolution due to absorption of light and 
%   diffusion of heat through the cuboid based on output of runMonteCarlo.m.
%   Also calculates Arrhenius-based thermal damage.
%
%	Pay attention to the sections with headers that say "USER SPECIFIED:"
%	In those sections, you must fill in the parameters relevant for your simulation.
%
%   Input
%       name
%           the basename of the files saved by defineGeometry.m and runMonteCarlo.m
%
%   Displays
%       Geometry cuboid
%       3D temperature evolution during illumination and diffusion
%       Temperature evolution for the individual temperature sensors (x-y-plot)
%       3D temperature distributions after illumination and after diffusion
%       Geometry cuboid showing any thermal damage
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
%           movie file showing the temperature evolution. The geometry cuboid
%           is shown in the beginning of the video.
%
%   Requires
%       deleteDataFiles.m
%       calcdtmax.m
%       plotVolumetric.m
%       updateVolumetric.m
%       finiteElementHeatPropagator.mex (architecture specific)
%
%% Load data from defineGeometry.m and runMonteCarlo.m
load(['./Data/' name '.mat'],'G');
load(['./Data/' name '_MCoutput.mat'],'MCoutput');

%% USER SPECIFIED: Define simulation behavior
% Should the script run in silent mode? (disables overwrite prompt,
% command window text, progress indication and plot generation)
silentMode = false;

% Should MCmatlab use all available processors on Windows? Otherwise,
% one will be left unused. Useful for doing other work on the PC
% while simulations are running.
useAllCPUs = true;

% Should a movie of the temperature evolution be generated? Requires silentMode = false.
makemovie  = false;

%% USER SPECIFIED: Define parameters
P                = 10; % [W] Incident pulse peak power (in case of infinite plane waves, only the power incident upon the cuboid's top surface)
onduration       = 0.001; % [s] Pulse on-duration
offduration      = 0.004; % [s] Pulse off-duration
Temp             = 37*ones(size(G.M)); % [deg C] Initial temperature distribution
n_pulses         = 1; % Number of consecutive pulses, each with an illumination phase and a diffusion phase. If simulating only illumination or only diffusion, use n_pulses = 1.

plotTempLimits   = [37 105]; % [deg C], the expected range of temperatures, used only for setting the color scale in the plot
n_updates_on     = 30; % Number of times data is extracted for plots during each illumination phase . Must be at least 1.
n_updates_off    = 120; % Number of times data is extracted for plots during each diffusion phase. Must be at least 1.

%% Check for preexisting files
if(~silentMode && ~deleteDataFiles(name)); return; end

%% Determine remaining parameters
if n_pulses ~= 1 && (onduration == 0 || offduration == 0)
    n_pulses = 1; % If either pulse on-duration or off-duration is 0, it doesn't make sense to talk of multiple pulses
    fprintf('\nWarning: Number of pulses changed to 1 since either on-duration or off-duration is 0.\n\n');
end

if(silentMode)
    n_updates_on  = 1;
    n_updates_off = 1;
end

[nx,ny,nz] = size(G.M);
nM = length(G.mediaProperties); % Number of different media in simulation

HC = G.dx*G.dy*G.dz*[G.mediaProperties.VHC]; % Heat capacity array. Element i is the HC of medium i in the reduced mediaProperties list.
TC = [G.mediaProperties.TC]; % Thermal conductivity array. Element i is the TC of medium i in the reduced mediaProperties list.

%% Prepare Arrhenius thermal damage parameters
A = [G.mediaProperties.A];
E = [G.mediaProperties.E];
if(any(A)) % If non-zero Arrhenius data exists, prepare to calculate thermal damage.
    Omega = zeros(size(G.M));
else
    Omega = NaN;
end

%% Calculate time step
dtmax = calcdtmax(G.M,TC,HC,G.dx,G.dy,G.dz)/2;

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
dTperdeltaT = cat(3,dt/G.dx*G.dy*G.dz*TC_eff./(diag(HC)*ones(nM)),...
                    dt*G.dx/G.dy*G.dz*TC_eff./(diag(HC)*ones(nM)),...
                    dt*G.dx*G.dy/G.dz*TC_eff./(diag(HC)*ones(nM))); % Third dimension corresponds to the different directions of heat diffusion
clear TC_eff

%% Calculate temperature change due to absorbed heat per time step
mua_vec = [G.mediaProperties.mua];
dT_abs = mua_vec(G.M).*MCoutput.F*dt*G.dx*G.dy*G.dz*P./HC(G.M); % Temperature change from absorption per time step [deg C] (mua_vec(G.M).*F is a 3D matrix of normalized volumetric powers)
clear MCoutput

if(~silentMode)
    %% Prepare the temperature plot
    if(~ishandle(21))
        heatsimFigure = figure(21);
        heatsimFigure.Position = [40 80 1100 650];
    else
        heatsimFigure = figure(21);
    end
    clf;
    heatsimFigure.Name = 'Temperature evolution';
    plotVolumetric(G.x,G.y,G.z,Temp,'MCmatlab');
    h_title = title(['Temperature evolution, t = ' num2str(timeVector(1),'%#.2g') ' s']);
    caxis(plotTempLimits); % User-defined color scale limits

    %% Plot the geometry to allow the user to select temperature sensor locations
    if(~ishandle(22))
        geometryFigure = figure(22);
        geometryFigure.Position = [40 80 1100 650];
    else
        geometryFigure = figure(22);
    end
    clf;
    geometryFigure.Name = 'Geometry illustration';
    plotVolumetric(G.x,G.y,G.z,G.M,'MCmatlab_GeometryIllustration',G.mediaProperties);
    title('Geometry illustration');

    datacursormode on;
    dataCursorHandle = datacursormode(geometryFigure);

    fprintf('Place one or more (shift+click) temperature sensor in the geometry,\n and hit (almost) any key to continue.\nYou can also resize the temperature figure now and position the slices,\n which is useful for movie generation.\n');
    pause
    cursorInfo = getCursorInfo(dataCursorHandle);
    numTemperatureSensors = length(cursorInfo);

    if numTemperatureSensors
        for temperatureSensorIndex = numTemperatureSensors:-1:1
            dataCursorPosition = cursorInfo(temperatureSensorIndex).Position;
            [~, dCPz] = min(abs(G.z-dataCursorPosition(3)));
            [~, dCPy] = min(abs(G.y-dataCursorPosition(2)));
            [~, dCPx] = min(abs(G.x-dataCursorPosition(1)));
            temperatureSensorPosition(temperatureSensorIndex) = ...
                sub2ind(size(G.M),dCPx,dCPy,dCPz);
        end
        temperatureSensorPosition = fliplr(temperatureSensorPosition); % Reverse array so the order of the sensors is the same as the order they were defined in.
        temperatureSensor = NaN(numTemperatureSensors,length(timeVector));
        temperatureSensor(:,1) = Temp(temperatureSensorPosition);

        temperatureSensorMedia = {G.mediaProperties(G.M(temperatureSensorPosition)).name};

        if(~ishandle(23))
            temperatureSensorFigure = figure(23);
            temperatureSensorFigure.Position = [40 80 1100 650];
        else
            temperatureSensorFigure = figure(23);
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
        legend(temperatureSensorMedia);
    else
        temperatureSensor = [];
    end

    if numTemperatureSensors; figure(temperatureSensorFigure); end

    fprintf('[nx,ny,nz]=[%d,%d,%d]. Number of pulses is %d.\nIllumination on for %d steps and off for %d steps in each pulse. Step size is %0.2e s.\n',nx,ny,nz,n_pulses,nt_on,nt_off,dt);

    if(makemovie) % Make a temporary figure showing the geometry illustration to put into the beginning of the movie
        tempFigure = figure;
        plotVolumetric(G.x,G.y,G.z,G.M,'MCmatlab_GeometryIllustration',G.mediaProperties);
        title('Geometry illustration');
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
else
    temperatureSensor = [];
end

%% Put heatSim parameters into a struct
heatSimParameters = struct('M',G.M-1,'Adt',A*dt,'E',E,'dTperdeltaT',dTperdeltaT,'dT_abs',dT_abs,'useAllCPUs',useAllCPUs); % Contents of G.M have to be converted from Matlab's 1-based indexing to C's 0-based indexing.

%% Simulate heat transfer
if(~silentMode); tic; end
for j=1:n_pulses
    if(~silentMode)
        if nt_on
            fprintf(['Illuminating pulse #' num2str(j) '... \n' repmat('-',1,n_updates_on)]);
        else
            fprintf(['Diffusing heat... \n' repmat('-',1,n_updates_off)]);
        end
        drawnow;
    end
    for i = 2:length(nt_vec)
        heatSimParameters.lightsOn = (nt_vec(i) <= nt_on);
        heatSimParameters.steps = nt_vec(i)-nt_vec(i-1);
        frameidx = i+(j-1)*(length(nt_vec)-1);
        [Temp,Omega] = finiteElementHeatPropagator(Temp,Omega,heatSimParameters);
        
        if(nt_vec(i) == nt_on && n_pulses == 1)
            Temp_illum = Temp;
        end
        
        if(~silentMode)
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
                fprintf('\b\b Done\n');
                if nt_off
                    fprintf(['Diffusing heat... \n' repmat('-',1,n_updates_off)]);
                end
            elseif nt_vec(i) == nt_on + nt_off
                fprintf('\b\b Done\n');
            end
        end
    end
end
if(~silentMode); toc; end
clear finiteElementHeatPropagator; % Unload finiteElementHeatPropagator MEX file so it can be modified externally again

%% Plot and save results
save(['./Data/' name '_heatSimoutput.mat'],'temperatureSensor','timeVector','P','onduration','offduration','Temp','Omega','n_pulses');

if(nt_on && nt_off && n_pulses == 1)
    save(['./Data/' name '_heatSimoutput.mat'],'Temp_illum','-append');
end

if(~silentMode)
    figure(heatsimFigure);
    if ~nt_on
        heatsimFigure.Name = 'Temperature after diffusion';
        title('Temperature after diffusion');
    elseif ~nt_off
        heatsimFigure.Name = 'Temperature after illumination';
        title('Temperature after illumination');
    elseif n_pulses == 1
        heatsimFigure.Name = 'Temperature after diffusion';
        title('Temperature after diffusion')
        if(~ishandle(24))
            illumFigure = figure(24);
            illumFigure.Position = [40 80 1100 650];
        else
            illumFigure = figure(24);
        end
        clf;
        illumFigure.Name = 'Temperature after illumination';
        plotVolumetric(G.x,G.y,G.z,Temp_illum,'MCmatlab');
        title('Temperature after illumination');
        caxis(plotTempLimits); % User-defined color scale limits
    else
        title(['Temperature after ' num2str(n_pulses) ' pulses']);
    end

    if ~isnan(Omega(1))
        if(~ishandle(25))
            damageFigure = figure(25);
            damageFigure.Position = [40 80 1100 650];
        else
            damageFigure = figure(25);
        end
        damageFigure.Name = 'Thermal damage illustration';
        M_damage = G.M;
        M_damage(Omega > 1) = nM + 1;
        G.mediaProperties(nM + 1).name = 'damage';
        plotVolumetric(G.x,G.y,G.z,M_damage,'MCmatlab_GeometryIllustration',G.mediaProperties);
        title('Thermal damage illustration');
        fprintf('%.2e cm^3 was thermally damaged.\n',G.dx*G.dy*G.dz*sum(sum(sum(Omega > 1))));
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
end
