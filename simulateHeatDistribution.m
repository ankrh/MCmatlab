function simulateHeatDistribution(name)
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
%       Geometry cuboid showing positions of specified temperature sensors
%       3D temperature evolution during illumination and diffusion
%       Temperature evolution for the individual temperature sensors (x-y-plot)
%       Geometry cuboid showing any thermal damage
%
%   Output
%       ./Data/[name]_heatSimOutput.mat
%           file containing the input values and the temperature sensor values at
%           each timestep
%
%       ./Data/[name]_heatSimOutput.mp4
%           movie file showing the temperature evolution. The geometry cuboid
%           is shown in the beginning of the video.
%
%   Requires
%       deleteDataFiles.m
%       calcdtmax.m
%       convertTempSensorPos.m
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
makemovie  = true;

% Would you like the cuboid boundaries to be insulated or heat-sinked?
% 0: Insulating boundaries
% 1: Constant-temperature boundaries (heat-sinked)
heatBoundaryType = 0;

%% USER SPECIFIED: Define parameters
P                = 4; % [W] Incident pulse peak power (in case of infinite plane waves, only the power incident upon the cuboid's top surface)
durationOn       = 0.001; % [s] Pulse on-duration
durationOff      = 0.004; % [s] Pulse off-duration
Temp             = 37*ones(size(G.M)); % [deg C] Initial temperature distribution
nPulses          = 5; % Number of consecutive pulses, each with an illumination phase and a diffusion phase. If simulating only illumination or only diffusion, use n_pulses = 1.

plotTempLimits   = [37 105]; % [deg C], the expected range of temperatures, used only for setting the color scale in the plot
nUpdates         = 100; % Number of times data is extracted for plots during each pulse. A minimum of 1 update is performed in each phase (2 for each pulse consisting of an illumination phase and a diffusion phase)

% Set the starting relative slice positions [x y z] for the 3D plots on a
% scale from 0 to 1. Especially relevant for movie generation. As an
% example, [0 1 0.5] puts slices at the lowest x value, the highest y value
% and the halfway z value.
slicePositions   = [.5 0.6 1];

% Programmatical temperature sensor positions, specified as a matrix where
% each row shows a temperature sensor's absolute [x y z] coordinates. Leave
% the matrix empty ([]) to disable temperature sensors.
tempSensorPositions = [0 0 0.038
                       0 0 0.04
                       0 0 0.042
                       0 0 0.044];

%% Check for preexisting files
if(~silentMode && ~deleteDataFiles(name)); return; end

%% Determine remaining parameters
if nPulses ~= 1 && (durationOn == 0 || durationOff == 0)
    nPulses = 1; % If either pulse on-duration or off-duration is 0, it doesn't make sense to talk of multiple pulses
    fprintf('\nWarning: Number of pulses changed to 1 since either on-duration or off-duration is 0.\n\n');
end

if(silentMode)
    nUpdatesOn  = 1;
    nUpdatesOff = 1;
end

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

if durationOn ~= 0
    nUpdatesOn = max(1,round(nUpdates*durationOn/(durationOn+durationOff)));
    nTsPerUpdateOn = ceil(durationOn/nUpdatesOn/dtmax); % Number of time steps per update in illumination phase
    dtOn = durationOn/nUpdatesOn/nTsPerUpdateOn; % Time step size during illumination
else
    nUpdatesOn = 0;
    nTsPerUpdateOn = 1;
    dtOn = 1;
end

if durationOff ~= 0
    nUpdatesOff = max(1,nUpdates - nUpdatesOn);
    nTsPerUpdateOff = ceil(durationOff/nUpdatesOff/dtmax); % Number of time steps per update in diffusion phase
    dtOff = durationOff/nUpdatesOff/nTsPerUpdateOff; % Time step size during diffusion
else
    nUpdatesOff = 0;
    nTsPerUpdateOff = 1;
    dtOff = 1;
end

updatesTimeVector = [0 , ((durationOn+durationOff)*repelem(0:(nPulses-1),nUpdatesOn+nUpdatesOff) + repmat([dtOn*nTsPerUpdateOn*(1:nUpdatesOn) , (durationOn + (dtOff*nTsPerUpdateOff*(1:nUpdatesOff)))],1,nPulses))];

sensorsTimeVector = [0 , ((durationOn+durationOff)*repelem(0:(nPulses-1),nUpdatesOn*nTsPerUpdateOn+nUpdatesOff*nTsPerUpdateOff) + repmat([dtOn*(1:nUpdatesOn*nTsPerUpdateOn) , (durationOn + (dtOff*(1:nUpdatesOff*nTsPerUpdateOff)))],1,nPulses))];

%% Calculate proportionality between voxel-to-voxel temperature difference DeltaT and time step temperature change dT
TC_eff = 2*(TC'*TC)./(ones(length(TC))*diag(TC)+diag(TC)*ones(length(TC))); % Same as TC_eff(i,j) = 2*TC_red(i)*TC_red(j)/(TC_red(i)+TC_red(j)) but without for loops
TC_eff(isnan(TC_eff)) = 0; % Neighboring insulating voxels return NaN but should just be 0

dTdtperdeltaT = cat(3,TC_eff./(diag(HC)*ones(nM))/G.dx*G.dy*G.dz,...
                      TC_eff./(diag(HC)*ones(nM))*G.dx/G.dy*G.dz,...
                      TC_eff./(diag(HC)*ones(nM))*G.dx*G.dy/G.dz); % [deg C /(s*deg C)] Time derivative of voxel temperature per voxel-to-voxel temperature difference. Third dimension corresponds to the different directions of heat diffusion (x, y and z)
clear TC_eff

%% Calculate temperature change due to absorbed heat per time step
mua_vec = [G.mediaProperties.mua];
dTdt_abs = mua_vec(G.M).*MCoutput.F*G.dx*G.dy*G.dz*P./HC(G.M); % [deg C/s] Time derivative of voxel temperature due to absorption. (mua_vec(G.M).*F is a 3D matrix of normalized volumetric powers)
clear MCoutput

%% Convert the temperature sensor positions to interpolation corner indices and weights used in the MEX function
numTemperatureSensors = size(tempSensorPositions,1);
if(numTemperatureSensors)
    if(any(abs(tempSensorPositions(:,1))               > G.dx*G.nx/2) || ...
       any(abs(tempSensorPositions(:,2))               > G.dy*G.ny/2) || ...
       any(abs(tempSensorPositions(:,3)) - G.dz*G.nz/2 > G.dz*G.nz/2))
        error('Error: A temperature sensor is outside the cuboid.');
    end

    [tempSensorCornerIdxs , tempSensorInterpWeights] = convertTempSensorPos(tempSensorPositions,G.x,G.y,G.z);
else
    tempSensorCornerIdxs = [];
    tempSensorInterpWeights = [];
end
sensorTemps = [];

if(~silentMode)
    if(numTemperatureSensors)
        %% Plot the geometry showing the temperature sensor locations
        if(~ishandle(22))
            geometryFigure = figure(22);
            geometryFigure.Position = [40 80 1100 650];
        else
            geometryFigure = figure(22);
        end
        geometryFigure.Name = 'Temperature sensor illustration';
        plotVolumetric(G.x,G.y,G.z,G.M,'MCmatlab_GeometryIllustration',G.mediaProperties,'slicePositions',slicePositions);
        title('Temperature sensor illustration');

        for i=numTemperatureSensors:-1:1
            sensorLabels{i,1} = num2str(i);
        end
        text(tempSensorPositions(:,1),tempSensorPositions(:,2),tempSensorPositions(:,3),sensorLabels,'HorizontalAlignment','center','FontSize',18)
        drawnow;
    end

    if durationOn ~= 0
        fprintf('Illumination phase consists of %d steps of %0.2e s.\n',nUpdatesOn*nTsPerUpdateOn,dtOn);
    end
    if durationOff ~= 0
        fprintf('Diffusion phase consists of %d steps of %0.2e s.\n',nUpdatesOff*nTsPerUpdateOff,dtOff);
    end
    
    if(~ishandle(21))
        heatsimFigure = figure(21);
        heatsimFigure.Position = [40 80 1100 650];
    else
        heatsimFigure = figure(21);
    end
    
    if(makemovie) % Make a temporary figure showing the geometry illustration to put into the beginning of the movie
        plotVolumetric(G.x,G.y,G.z,G.M,'MCmatlab_GeometryIllustration',G.mediaProperties,'slicePositions',slicePositions);
        title('Geometry illustration');
        drawnow;
        movieframes(1) = getframe(heatsimFigure);
    end
    %% Prepare the temperature plot
    heatsimFigure.Name = 'Temperature evolution';
    plotVolumetric(G.x,G.y,G.z,Temp,'MCmatlab','slicePositions',slicePositions);
    h_title = title('Temperature evolution, t = 0 s');
    caxis(plotTempLimits); % User-defined color scale limits
    if(makemovie); movieframes(2) = getframe(heatsimFigure); end
end

%% Put heatSim parameters into a struct
heatSimParameters = struct('M',G.M-1,'A',A,'E',E,'dTdtperdeltaT',dTdtperdeltaT,'dTdt_abs',dTdt_abs,'useAllCPUs',useAllCPUs,'heatBoundaryType',heatBoundaryType,'tempSensorCornerIdxs',tempSensorCornerIdxs,'tempSensorInterpWeights',tempSensorInterpWeights); % Contents of G.M have to be converted from Matlab's 1-based indexing to C's 0-based indexing.
clear dTdtperdeltaT dTdt_abs % These 3D matrices are large and no longer needed since they have been copied to the struct, so they are cleared to conserve memory

%% Simulate heat transfer
if(~silentMode); tic; end
for j=1:nPulses
    if durationOn ~= 0
        if(~silentMode)
            fprintf(['Illuminating pulse #' num2str(j) '... \n' repmat('-',1,nUpdatesOn)]);
            drawnow;
        end
        
        heatSimParameters.lightsOn = true;
        heatSimParameters.steps = nTsPerUpdateOn;
        heatSimParameters.dt = dtOn;
        for i = 1:nUpdatesOn
            [Temp,Omega,newSensorTemps] = finiteElementHeatPropagator(Temp,Omega,heatSimParameters);
            
            if isempty(sensorTemps)
                sensorTemps = newSensorTemps;
            else
                sensorTemps = [sensorTemps(:,1:end-1) newSensorTemps];
            end
            
            if(~silentMode)
                fprintf(1,'\b');
                updateIdx = i+(j-1)*(nUpdatesOn+nUpdatesOff); % Index of the update
                updateVolumetric(heatsimFigure,Temp);
                h_title.String = ['Temperature evolution, t = ' num2str(updatesTimeVector(updateIdx+1),'%#.2g') ' s'];
                drawnow;
                
                if(makemovie)
                    movieframes(updateIdx+2) = getframe(heatsimFigure);
                end
            end
        end
        
        if(~silentMode); fprintf('\b\b Done\n'); end
    end
    
    if durationOff ~= 0
        if(~silentMode)
            fprintf(['Diffusing heat... \n' repmat('-',1,nUpdatesOff)]);
            drawnow;
        end
        
        heatSimParameters.lightsOn = false;
        heatSimParameters.steps = nTsPerUpdateOff;
        heatSimParameters.dt = dtOff;
        for i = 1:nUpdatesOff
            [Temp,Omega,newSensorTemps] = finiteElementHeatPropagator(Temp,Omega,heatSimParameters);
            
            if isempty(sensorTemps)
                sensorTemps = newSensorTemps;
            else
                sensorTemps = [sensorTemps(:,1:end-1) newSensorTemps];
            end
            
            if(~silentMode)
                fprintf(1,'\b');
                updateIdx = nUpdatesOn+i+(j-1)*(nUpdatesOn+nUpdatesOff); % Index of the update
                updateVolumetric(heatsimFigure,Temp);
                h_title.String = ['Temperature evolution, t = ' num2str(updatesTimeVector(updateIdx+1),'%#.2g') ' s'];
                drawnow;
                
                if(makemovie)
                    movieframes(updateIdx+2) = getframe(heatsimFigure);
                end
            end
        end
        if(~silentMode); fprintf('\b\b Done\n'); end
    end
end
if(~silentMode); toc; end
clear finiteElementHeatPropagator; % Unload finiteElementHeatPropagator MEX file so it can be modified externally again
clear Temp heatSimParameters;

%% Plot and save results
save(['./Data/' name '_heatSimoutput.mat'],'sensorsTimeVector','sensorTemps','P','durationOn','durationOff','Omega','nPulses');

if(~silentMode)
    if(numTemperatureSensors)
        if(~ishandle(23))
            temperatureSensorFigure = figure(23);
            temperatureSensorFigure.Position = [40 80 1100 650];
        else
            temperatureSensorFigure = figure(23);
        end
        clf;
        temperatureSensorFigure.Name = 'Temperature sensors';
        plot(sensorsTimeVector,sensorTemps,'LineWidth',2);
        set(gca,'FontSize',16);
        xlabel('Time [sec]')
        ylabel('Temperature [deg C]')
        title('Temperature sensors')
        xlim(sensorsTimeVector([1 end]));
        legend(sensorLabels,'Location','best');
        grid on;grid minor;
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
        plotVolumetric(G.x,G.y,G.z,M_damage,'MCmatlab_GeometryIllustration',G.mediaProperties,'slicePositions',slicePositions);
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
