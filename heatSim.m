function [temperatureSensor, timeVector] = heatSim()
%% User defined parameters

directoryPath = 'exec/';
myname = 'dentin_sim_850';
Winc            = 1; % [W] Incident pulse peak power (in case of infinite plane waves, only the power incident upon the box' top area)
onduration      = 1; % [s] Pulse on-duration
offduration     = 1; % [s] Pulse off-duration
initialTemp     = 0; % [deg C]
plotTempLimits  = [0 1e-3]; % [deg C], the expected range of temperatures, used only for setting the color scale in the plot
extremeprecision = false; % false: normal time step resolution (dt = dtmax/2), true: high time step resolution (dt = dtmax/100)
n_updates       = 100; % Number of times data is extracted for plots during illumination phase and diffusion phase. Must be at least 1.

%% Load data
load([directoryPath myname '_HS.mat']);

%% Parameters that are set automatically
[nx,ny,nz] = size(T);

x  = ((0:nx-1)-(nx-1)/2)*dx;
y  = ((0:ny-1)-(ny-1)/2)*dy;
z  = ((0:nz-1)+1/2)*dz;

%% Define reduced matrices and vectors
tissueMap = zeros(1,length(tissueList),'uint8');
tissueMap(unique(T)) = 1:length(unique(T));
tissueList = tissueList(unique(T));
T = tissueMap(T); % Reduced tissue matrix, using only numbers from 1 up to the number of tissues actually used
HC = dx*dy*dz*[tissueList.VHC]; % Reduced heat capacity vector. Element i is the HC of tissue i in the reduced tissue list.
TC = [tissueList.TC]; % Reduced thermal conductivity vector. Element i is the TC of tissue i in the reduced tissue list.

%% Calculate time step
individual_dtmax = calcdtmax(T,TC,HC,dx,dy,dz);
if extremeprecision
    dtmax = min(individual_dtmax(:))/100;
else
    dtmax = min(individual_dtmax(:))/2;
end
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
timeVector = nt_vec*dt; % For temperature sensor plotting

%% Calculate proportionality between temperature difference DeltaT and diffused heat per time step dQ
TC_eff = 2*(TC'*TC)./(ones(length(TC))*diag(TC)+diag(TC)*ones(length(TC))); % Same as TC_eff(i,j) = 2*TC_red(i)*TC_red(j)/(TC_red(i)+TC_red(j)) but without for loops
TC_eff(isnan(TC_eff)) = 0; % Neighboring insulating voxels return NaN but should just be 0
dQperdeltaT = cat(3,dt/dx*dy*dz*TC_eff,dt*dx/dy*dz*TC_eff,dt*dx*dy/dz*TC_eff); % Third dimension corresponds to the different directions of heat diffusion
clear TC_eff

%% Calculate absorbed heat per time step
dQ_abs = NVP*dt*dx*dy*dz*Winc; % Heat from absorption per time step [J]
clear NVP

%% Initialize temperature
Temp = initialTemp*ones(size(T));

%% Plots to visualize tissue properties

% figure(1);clf;
% plotVolumetric(x,y,z,HC(Tissue)/(dx*dy*dz));
% title('Volumetric heat capacity [J/cm^3/K]');
% figure(2);clf;
% plotVolumetric(x,y,z,TC(Tissue));
% title('Thermal conductivity [W/cm/K]');

tissueFigure = figure(3); clf;
plotVolumetric(x,y,z,T,tissueList);

datacursormode on;
dataCursorHandle = datacursormode(tissueFigure);

fprintf('Place one or more (shift+click) temperature sensor in the tissue,\n and hit (almost) any key to continue...\n')
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
    temperatureSensor = NaN(numTemperatureSensors,length(nt_vec));
    temperatureSensor(:,1) = Temp(temperatureSensorPosition);
    
    temperatureSensorFigure = figure(4); clf;
    temperaturePlot = axes('XLim',timeVector([1, end]));
    xlabel('Time [sec]')
    ylabel('Temperature [?C]')
    title('Temperature evolution')
    hold on
end

%% Prepare the temperature plot
heatsimFigure = figure(5);
plotVolumetric(x,y,z,Temp);
title('Temperature evolution');
caxis(plotTempLimits); % User-defined color scale limits

fprintf('Adjust the visualisation to your needs, and press (almost) any key to start simulation...\n')
pause

%% Heat Transfer Simulation during illumination
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
    Temp = propagator(nt_vec(i)-nt_vec(i-1),Temp,T-1,dQperdeltaT,(nt_vec(i) <= nt_on)*dQ_abs,HC); % Arguments (nt,[[[Temp]]],[[[T]]],[[[dQperdeltaT]]],[[[dQ_abs]]],[HC]). Contents of T has to be converted from Matlab's 1-based indexing to C's 0-bsed indexing.
    
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

figure(heatsimFigure)
if ~nt_on
    title('Temperature after diffusion');
elseif ~nt_off
    title('Temperature after illumination');
else
    title('Temperature after diffusion')
    figure(6);
    plotVolumetric(x,y,z,Temp_illum);
    title('Temperature after illumination');
end
end


