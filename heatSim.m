function [temperatureSensor, timeVector] = heatSim(name)
%% Load data from makeTissue.m and MonteCarlo.m
load(['./Data/' name '.mat']);
load(['./Data/' name '_MCoutput.mat'],'F');

%% Define parameters (user-specified)
Winc             = 1; % [W] Incident pulse peak power (in case of infinite plane waves, only the power incident upon the box' top area)
onduration       = 0.5; % [s] Pulse on-duration
offduration      = 0.5; % [s] Pulse off-duration
Temp             = zeros(size(T)); % [deg C] Initial temperature distrubution

plotTempLimits   = [0 16]; % [deg C], the expected range of temperatures, used only for setting the color scale in the plot
extremeprecision = false; % false: normal time step resolution (dt = dtmax/2), true: high time step resolution (dt = dtmax/100)
n_updates        = 20; % Number of times data is extracted for plots during illumination phase and diffusion phase. Must be at least 1.

%% Determine remaining parameters
[nx,ny,nz] = size(T);
nT = length(tissueList); % Number of different tissues in simulation

dx = x(2) - x(1);
dy = y(2) - y(1);
dz = z(2) - z(1);

HC = dx*dy*dz*[tissueList.VHC]; % Heat capacity array. Element i is the HC of tissue i in the reduced tissue list.
TC = [tissueList.TC]; % Thermal conductivity array. Element i is the TC of tissue i in the reduced tissue list.

%% Calculate time step
if extremeprecision
    dtmax = calcdtmax(T,TC,HC,dx,dy,dz)/100;
else
    dtmax = calcdtmax(T,TC,HC,dx,dy,dz)/2;
end

% figure(7);
% plotVolumetric(x,y,z,individual_dtmax); % Optional plot to visualize each individual voxel's required time step limitation
% title('Maximum time step allowed for each voxel [s]')

clear individual_dtmax

if onduration ~= 0
    nt_on = max(n_updates,ceil(onduration/dtmax)); % Number of time steps with illumination
    dt = onduration/nt_on; % Time step size
    if offduration ~= 0
        nt_off = max(n_updates,ceil(offduration/dt)); % Number of time steps without illumination
    else
        nt_off = 0;
    end
else
    nt_on = 0;
    nt_off = max(n_updates,ceil(offduration/dtmax)); % Number of time steps without illumination
    dt = offduration/nt_off; % Time step size
end

nt_vec = unique(floor([linspace(0,nt_on,n_updates+1) linspace(nt_on + nt_off/n_updates,nt_on + nt_off,n_updates)])); % Vector of nt's on which the plots should be updated

%% Calculate proportionality between voxel-to-voxel temperature difference DeltaT and time step temperature change dT
TC_eff = 2*(TC'*TC)./(ones(length(TC))*diag(TC)+diag(TC)*ones(length(TC))); % Same as TC_eff(i,j) = 2*TC_red(i)*TC_red(j)/(TC_red(i)+TC_red(j)) but without for loops
TC_eff(isnan(TC_eff)) = 0; % Neighboring insulating voxels return NaN but should just be 0
dTperdeltaT = cat(3,dt/dx*dy*dz*TC_eff./(diag(HC)*ones(nT)),dt*dx/dy*dz*TC_eff./(diag(HC)*ones(nT)),dt*dx*dy/dz*TC_eff./(diag(HC)*ones(nT))); % Third dimension corresponds to the different directions of heat diffusion
clear TC_eff

%% Calculate temperature change due to absorbed heat per time step
mua_vec = [tissueList.mua];
dT_abs = mua_vec(T).*F*dt*dx*dy*dz*Winc./HC(T); % Temperature change from absorption per time step [deg C] (mua_vec(T).*F is a 3D matrix of normalized volumetric powers)
clear F

%% Make plots to visualize tissue properties

% figure(2);clf;
% plotVolumetric(x,y,z,HC(T)/(dx*dy*dz));
% title('Volumetric heat capacity [J/cm^3/K]');
% figure(3);clf;
% plotVolumetric(x,y,z,TC(T));
% title('Thermal conductivity [W/cm/K]');

tissueFigure = figure(1); clf;
plotVolumetric(x,y,z,T,tissueList);
title('Tissue type illustration');

datacursormode on;
dataCursorHandle = datacursormode(tissueFigure);

fprintf('Place one or more (shift+click) temperature sensor in the tissue,\n and hit (almost) any key to continue...\n')
pause
cursorInfo = getCursorInfo(dataCursorHandle);
numTemperatureSensors = length(cursorInfo);

if numTemperatureSensors
    timeVector = nt_vec*dt; % For temperature sensor plotting
    for temperatureSensorIndex = numTemperatureSensors:-1:1
        dataCursorPosition = cursorInfo(temperatureSensorIndex).Position;
        [~, dCPz] = min(abs(z-dataCursorPosition(3)));
        [~, dCPy] = min(abs(y-dataCursorPosition(2)));
        [~, dCPx] = min(abs(x-dataCursorPosition(1)));
        temperatureSensorPosition(temperatureSensorIndex) = ...
            sub2ind(size(T),dCPx,dCPy,dCPz);
    end
    temperatureSensor = NaN(numTemperatureSensors,length(nt_vec));
    temperatureSensor(:,1) = Temp(temperatureSensorPosition);
    
    temperatureSensorTissues = {tissueList(T(temperatureSensorPosition)).name};
    
    temperatureSensorFigure = figure(4); clf;
    temperaturePlot = axes('XLim',timeVector([1, end]));
    xlabel('Time [sec]')
    ylabel('Temperature [?C]')
    title('Temperature evolution')
    hold on
    plot(temperaturePlot,timeVector,temperatureSensor);
    legend(temperatureSensorTissues);
else
    temperatureSensor = [];
    timeVector = [];
end

%% Prepare the temperature plot
heatsimFigure = figure(5);
plotVolumetric(x,y,z,Temp);
title('Temperature evolution');
caxis(plotTempLimits); % User-defined color scale limits

fprintf('Adjust the visualisation to your needs, and press (almost) any key to start simulation...\n')
pause

%% Simulate heat transfer
figure(heatsimFigure)
if numTemperatureSensors; figure(temperatureSensorFigure); end

fprintf('[nx,ny,nz]=[%d,%d,%d]\nIllumination on for %d steps and off for %d steps. Step size is %0.2e s.\n',nx,ny,nz,nt_on,nt_off,dt);

tic
if nt_on
    fprintf(['Illuminating... \n' repmat('-',1,n_updates)]);
else
    fprintf(['Diffusing heat... \n' repmat('-',1,n_updates)]);
end
drawnow;
for i = 2:length(nt_vec)
    Temp = propagator(nt_vec(i)-nt_vec(i-1),Temp,T-1,dTperdeltaT,(nt_vec(i) <= nt_on)*dT_abs); % Arguments (nt,[[[Temp]]],[[[T]]],[[[dTperdeltaT]]],[[[dT_abs]]]). Contents of T have to be converted from Matlab's 1-based indexing to C's 0-based indexing.
    
    if numTemperatureSensors
        temperatureSensor(:,i) = Temp(temperatureSensorPosition);
    end
    
    fprintf(1,'\b');
    updateVolumetric(heatsimFigure,Temp);
    if numTemperatureSensors
        cla(temperaturePlot)
        plot(temperaturePlot,timeVector,temperatureSensor);
    end
    
    if nt_vec(i) == nt_on
        Temp_illum = Temp;
        fprintf('\b\b Done\n');
        if nt_off
            fprintf(['Diffusing heat... \n' repmat('-',1,n_updates)]);
        end
    elseif nt_vec(i) == nt_on + nt_off
        fprintf('\b\b Done\n');
    end
end
toc;

%% Plot and save results
figure(heatsimFigure)
save(['./Data/' name '_heatSimoutput.mat'],'temperatureSensor','timeVector','Winc','onduration','offduration','Temp');
if ~nt_on
    title('Temperature after diffusion');
elseif ~nt_off
    title('Temperature after illumination');
else
    title('Temperature after diffusion')
    figure(6);
    plotVolumetric(x,y,z,Temp_illum);
    title('Temperature after illumination');
    save(['./Data/' name '_heatSimoutput.mat'],'Temp_illum','-append');
end
fprintf('./Data/%s_heatSimoutput.mat saved\n',name);
end