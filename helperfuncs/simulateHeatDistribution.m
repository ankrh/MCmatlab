function HSoutput = simulateHeatDistribution(HSinput)
%   Simulates temperature evolution due to absorption of light and 
%   diffusion of heat through the cuboid based on output of runMonteCarlo.m.
%   Also calculates Arrhenius-based thermal damage.
%
%   Output
%       ./Data/[name]_heatSimOutput.mp4
%           movie file showing the temperature evolution. The geometry cuboid
%           is shown in the beginning of the video.
%
%   Requires
%       calcdtmax.m
%       convertTempSensorPos.m
%       plotVolumetric.m
%       updateVolumetric.m
%       finiteElementHeatPropagator.mex (architecture specific)
%
%	See also runMonteCarlo, plotMCmatlabHeat

%%%%%
%   Copyright 2017, 2018 by Dominik Marti and Anders K. Hansen, DTU Fotonik
%   Work on this heat solver was started by Rasmus L. Pedersen & Mathias Christensen, DTU Fotonik
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

G = HSinput.G;

Temp = HSinput.initialTemp*ones(size(G.M)); % Initial temperature distribution
nM = length(G.mediaProperties); % Number of different media in simulation

HC = G.dx*G.dy*G.dz*[G.mediaProperties.VHC]; % Heat capacity array. Element i is the HC of medium i in the reduced mediaProperties list.
TC = [G.mediaProperties.TC]; % Thermal conductivity array. Element i is the TC of medium i in the reduced mediaProperties list.

%% Prepare Arrhenius thermal damage parameters
A = [G.mediaProperties.A];
E = [G.mediaProperties.E];
if(any(A)) % If non-zero Arrhenius data exists, prepare to calculate thermal damage.
    HSoutput.Omega = zeros(size(G.M));
else
    HSoutput.Omega = NaN;
end

%% Calculate amount and duration of time steps in each phase
% If either pulse on-duration or off-duration is 0, it doesn't make sense to talk of multiple pulses
if HSinput.nPulses ~= 1 && (HSinput.durationOn == 0 || HSinput.durationOff == 0)
    HSinput.nPulses = 1;
    fprintf('\nWarning: Number of pulses changed to 1 since either on-duration or off-duration is 0.\n\n');
end

if(HSinput.silentMode)
    HSinput.nUpdates = 1;
end

dtmax = calcdtmax(G.M,TC,HC,G.dx,G.dy,G.dz)/2; % Highest allowable time step duration

if HSinput.durationOn ~= 0
    nUpdatesOn = max(1,round(HSinput.nUpdates*HSinput.durationOn/(HSinput.durationOn+HSinput.durationOff)));
    nTsPerUpdateOn = ceil(HSinput.durationOn/nUpdatesOn/dtmax); % Number of time steps per update in illumination phase
    dtOn = HSinput.durationOn/nUpdatesOn/nTsPerUpdateOn; % Time step size during illumination
else
    nUpdatesOn = 0;
    nTsPerUpdateOn = 1;
    dtOn = 1;
end

if HSinput.durationOff ~= 0
    nUpdatesOff = max(1,HSinput.nUpdates - nUpdatesOn);
    nTsPerUpdateOff = ceil(HSinput.durationOff/nUpdatesOff/dtmax); % Number of time steps per update in diffusion phase
    dtOff = HSinput.durationOff/nUpdatesOff/nTsPerUpdateOff; % Time step size during diffusion
else
    nUpdatesOff = 0;
    nTsPerUpdateOff = 1;
    dtOff = 1;
end

if HSinput.durationEnd ~= 0
    nUpdatesEnd = max(1,round(HSinput.nUpdates*HSinput.durationEnd/(HSinput.durationOn+HSinput.durationOff)));
    nTsPerUpdateEnd = ceil(HSinput.durationEnd/nUpdatesEnd/dtmax); % Number of time steps per update in end relaxation phase
    dtEnd = HSinput.durationEnd/nUpdatesEnd/nTsPerUpdateEnd; % Time step size during end relaxation phase
else
    nUpdatesEnd = 0;
    nTsPerUpdateEnd = 1;
    dtEnd = 1;
end

% Array of times associated with the updates
updatesTimeVector = [0 , ((HSinput.durationOn+HSinput.durationOff)*repelem(0:(HSinput.nPulses-1),nUpdatesOn                + nUpdatesOff                ) + repmat([(1:nUpdatesOn)*nTsPerUpdateOn*dtOn , (HSinput.durationOn + ((1:nUpdatesOff)*nTsPerUpdateOff*dtOff))],1,HSinput.nPulses)) , (HSinput.durationOn + HSinput.durationOff)*HSinput.nPulses + (1:nUpdatesEnd)*nTsPerUpdateEnd*dtEnd];

% Array of times associated with the steps
HSoutput.sensorsTimeVector = [0 , ((HSinput.durationOn+HSinput.durationOff)*repelem(0:(HSinput.nPulses-1),nUpdatesOn*nTsPerUpdateOn + nUpdatesOff*nTsPerUpdateOff) + repmat([(1:nUpdatesOn*nTsPerUpdateOn)*dtOn , (HSinput.durationOn + ((1:nUpdatesOff*nTsPerUpdateOff)*dtOff))],1,HSinput.nPulses)) , (HSinput.durationOn + HSinput.durationOff)*HSinput.nPulses + (1:nUpdatesEnd*nTsPerUpdateEnd)*dtEnd];

%% Calculate proportionality between voxel-to-voxel temperature difference DeltaT and time step temperature change dT
TC_eff = 2*(TC'*TC)./(ones(length(TC))*diag(TC)+diag(TC)*ones(length(TC))); % Same as TC_eff(i,j) = 2*TC_red(i)*TC_red(j)/(TC_red(i)+TC_red(j)) but without for loops
TC_eff(isnan(TC_eff)) = 0; % Neighboring insulating voxels return NaN but should just be 0

dTdtperdeltaT = cat(3,TC_eff./(diag(HC)*ones(nM))/G.dx*G.dy*G.dz,...
                      TC_eff./(diag(HC)*ones(nM))*G.dx/G.dy*G.dz,...
                      TC_eff./(diag(HC)*ones(nM))*G.dx*G.dy/G.dz); % [deg C /(s*deg C)] Time derivative of voxel temperature per voxel-to-voxel temperature difference. Third dimension corresponds to the different directions of heat diffusion (x, y and z)
clear TC_eff

%% Calculate temperature change due to absorbed heat per time step
mua_vec = [G.mediaProperties.mua];
dTdt_abs = mua_vec(G.M).*HSinput.MCoutput.F*G.dx*G.dy*G.dz*HSinput.P./HC(G.M); % [deg C/s] Time derivative of voxel temperature due to absorption. (mua_vec(G.M).*F is a 3D matrix of normalized volumetric powers)

%% Convert the temperature sensor positions to interpolation corner indices and weights used in the MEX function
numTemperatureSensors = size(HSinput.tempSensorPositions,1);
if(numTemperatureSensors)
    if(any(abs(HSinput.tempSensorPositions(:,1))               > G.dx*G.nx/2) || ...
       any(abs(HSinput.tempSensorPositions(:,2))               > G.dy*G.ny/2) || ...
       any(abs(HSinput.tempSensorPositions(:,3)) - G.dz*G.nz/2 > G.dz*G.nz/2))
        error('Error: A temperature sensor is outside the cuboid.');
    end

    [tempSensorCornerIdxs , tempSensorInterpWeights] = convertTempSensorPos(HSinput.tempSensorPositions,G.x,G.y,G.z);
else
    tempSensorCornerIdxs = [];
    tempSensorInterpWeights = [];
end
HSoutput.sensorTemps = [];

%% Output some diagnostics info and prepare the temperature evolution plot. If making a movie, put a geometry illustration into the beginning of the movie.
if(~HSinput.silentMode)
    if HSinput.durationOn ~= 0
        fprintf('Illumination phase consists of %d steps of %0.2e s.\n',nUpdatesOn*nTsPerUpdateOn,dtOn);
    end
    if HSinput.durationOff ~= 0
        fprintf('Diffusion phase consists of %d steps of %0.2e s.\n',nUpdatesOff*nTsPerUpdateOff,dtOff);
    end
    if HSinput.durationEnd ~= 0
        fprintf('End phase consists of %d steps of %0.2e s.\n',nUpdatesEnd*nTsPerUpdateEnd,dtEnd);
    end
    
    if(HSinput.makemovie) % Make a temporary figure showing the geometry illustration to put into the beginning of the movie
        heatsimFigure = plotVolumetric(21,G.x,G.y,G.z,G.M,'MCmatlab_GeometryIllustration',G.mediaProperties,'slicePositions',HSinput.slicePositions);
        title('Geometry illustration');
        drawnow;
        movieframes(1) = getframe(heatsimFigure);
    end
    
    heatsimFigure = plotVolumetric(21,G.x,G.y,G.z,Temp,'MCmatlab_heat','slicePositions',HSinput.slicePositions);
    heatsimFigure.Name = 'Temperature evolution';
    h_title = title('Temperature evolution, t = 0 s');
    caxis(HSinput.plotTempLimits); % User-defined color scale limits
    if(HSinput.makemovie); movieframes(2) = getframe(heatsimFigure); end
end

%% Put heatSim parameters into a struct
heatSimParameters = struct('M',G.M-1,'A',A,'E',E,'dTdtperdeltaT',dTdtperdeltaT,'dTdt_abs',dTdt_abs,...
                'useAllCPUs',HSinput.useAllCPUs,'heatBoundaryType',HSinput.heatBoundaryType,...
                'tempSensorCornerIdxs',tempSensorCornerIdxs,'tempSensorInterpWeights',tempSensorInterpWeights); % Contents of G.M have to be converted from Matlab's 1-based indexing to C's 0-based indexing.

%% Simulate heat transfer
if(~HSinput.silentMode); tic; end
for j=1:HSinput.nPulses
	%% Illumination phase
    if HSinput.durationOn ~= 0
        if(~HSinput.silentMode)
            fprintf(['Illuminating pulse #' num2str(j) '... \n' repmat('-',1,nUpdatesOn)]);
            drawnow;
        end
        
        heatSimParameters.lightsOn = true;
        heatSimParameters.steps = nTsPerUpdateOn;
        heatSimParameters.dt = dtOn;
        for i = 1:nUpdatesOn
            [Temp,HSoutput.Omega,newSensorTemps] = finiteElementHeatPropagator(Temp,HSoutput.Omega,heatSimParameters);
            
            if isempty(HSoutput.sensorTemps)
                HSoutput.sensorTemps = newSensorTemps;
            else
                HSoutput.sensorTemps = [HSoutput.sensorTemps(:,1:end-1) newSensorTemps];
            end
            
            if(~HSinput.silentMode)
                fprintf(1,'\b');
                updateIdx = i+(j-1)*(nUpdatesOn+nUpdatesOff); % Index of the update
                updateVolumetric(heatsimFigure,Temp);
                h_title.String = ['Temperature evolution, t = ' num2str(updatesTimeVector(updateIdx+1),'%#.2g') ' s'];
                drawnow;
                
                if(HSinput.makemovie)
                    movieframes(updateIdx+2) = getframe(heatsimFigure);
                end
            end
        end
        
        if(~HSinput.silentMode); fprintf('\b\b Done\n'); end
    end
    
	%% Diffusion phase
    if HSinput.durationOff ~= 0
        if(~HSinput.silentMode)
            fprintf(['Diffusing heat... \n' repmat('-',1,nUpdatesOff)]);
            drawnow;
        end
        
        heatSimParameters.lightsOn = false;
        heatSimParameters.steps = nTsPerUpdateOff;
        heatSimParameters.dt = dtOff;
        for i = 1:nUpdatesOff
            [Temp,HSoutput.Omega,newSensorTemps] = finiteElementHeatPropagator(Temp,HSoutput.Omega,heatSimParameters);
            
            if isempty(HSoutput.sensorTemps)
                HSoutput.sensorTemps = newSensorTemps;
            else
                HSoutput.sensorTemps = [HSoutput.sensorTemps(:,1:end-1) newSensorTemps];
            end
            
            if(~HSinput.silentMode)
                fprintf(1,'\b');
                updateIdx = nUpdatesOn+i+(j-1)*(nUpdatesOn+nUpdatesOff); % Index of the update
                updateVolumetric(heatsimFigure,Temp);
                h_title.String = ['Temperature evolution, t = ' num2str(updatesTimeVector(updateIdx+1),'%#.2g') ' s'];
                drawnow;
                
                if(HSinput.makemovie)
                    movieframes(updateIdx+2) = getframe(heatsimFigure);
                end
            end
        end
        if(~HSinput.silentMode); fprintf('\b\b Done\n'); end
    end
end

%% End relaxation phase
if(~HSinput.silentMode)
	fprintf(['Diffusing heat in end relaxation phase... \n' repmat('-',1,nUpdatesEnd)]);
	drawnow;
end

heatSimParameters.lightsOn = false;
heatSimParameters.steps = nTsPerUpdateEnd;
heatSimParameters.dt = dtEnd;
for i = 1:nUpdatesEnd
	[Temp,HSoutput.Omega,newSensorTemps] = finiteElementHeatPropagator(Temp,HSoutput.Omega,heatSimParameters);

	if isempty(HSoutput.sensorTemps)
		HSoutput.sensorTemps = newSensorTemps;
	else
		HSoutput.sensorTemps = [HSoutput.sensorTemps(:,1:end-1) newSensorTemps];
	end

	if(~HSinput.silentMode)
		fprintf(1,'\b');
		updateIdx = i + HSinput.nPulses*(nUpdatesOn+nUpdatesOff); % Index of the update
		updateVolumetric(heatsimFigure,Temp);
		h_title.String = ['Temperature evolution, t = ' num2str(updatesTimeVector(updateIdx+1),'%#.2g') ' s'];
		drawnow;

		if(HSinput.makemovie)
			movieframes(updateIdx+2) = getframe(heatsimFigure);
		end
	end
end

if(~HSinput.silentMode)
    fprintf('\b\b Done\n');
    toc;
end

clear finiteElementHeatPropagator; % Unload finiteElementHeatPropagator MEX file so it can be modified externally again
clear Temp heatSimParameters;

%% Finalize and write movie
if(~HSinput.silentMode && HSinput.makemovie)
    movieframes = [repmat(movieframes(1),1,30) movieframes(1:end) repmat(movieframes(end),1,30)];
    writerObj = VideoWriter(['./Data/heatSimoutput.mp4'],'MPEG-4');
    writerObj.Quality = 100;
    open(writerObj);
    warning('off','MATLAB:audiovideo:VideoWriter:mp4FramePadded');
    writeVideo(writerObj,movieframes);
    warning('on','MATLAB:audiovideo:VideoWriter:mp4FramePadded');
    close(writerObj);
    fprintf('./Data/heatSimoutput.mp4 saved\n');
end
end
