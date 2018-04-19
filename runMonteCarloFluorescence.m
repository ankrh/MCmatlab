function runMonteCarloFluorescence(name)
% Script for simulating distribution and magnitude of fluorescence
% based on runMonteCarlo.m

% Created 2018-04-19 by Anders K. Hansen, DTU Fotonik

%% Load data from makeTissue.m and runMonteCarlo.m
load(['./Data/' name '.mat']);
load(['./Data/' name '_MCoutput.mat']);

MCinput_fluorescence = MCinput;

%% Define parameters (user-specified)
% Simulation duration
MCinput_fluorescence.simulationTime = 0.1;      % [min] time duration of the simulation

% Incident excitation power, relevant only if fluorescence saturations are finite
P = 2; % [W]

% Boundary type
% 0 = no boundaries
% 1 = escape at boundaries
% 2 = escape at surface only. No x, y, bottom z boundaries
MCinput_fluorescence.boundaryFlag = 1;

%% Set remaining parameters
MCinput_fluorescence.tissueList = tissueList_fluorescence;

%% Calculate 3D fluorescence source distribution, including saturation
mua_vec = [tissueList.mua]; % The tissues' excitation absorption coefficients
Y_vec = [tissueList.Y]; % The tissues' fluorescence power yields
sat_vec = [tissueList.sat]; % The tissues' fluorescence saturation fluence rates (intensity)
MCinput_fluorescence.sourceDistribution = Y_vec(T).*mua_vec(T)*P.*F./(1 + P*F./sat_vec(T)); % [W/cm^3]

%% Call Monte Carlo C script to get fluorescence distribution
I_fluorescence = mcxyz(MCinput_fluorescence); % This is an absolute fluence rate (intensity) quantity, unlike the "F" matrices which are actually fluence rates normalized to the incident power

%% Check for preexisting files
if(exist(['./Data/' name '_MCoutput_fluorescence.mat'],'file'))
    if(strcmp(questdlg('Computation results by this name already exist. Delete existing files?','Overwrite prompt','Yes','No, abort','Yes'),'No, abort'))
        fprintf('Aborted without saving data.\n');
        return;
    end
end

%% Save output and clear memory
save(['./Data/' name '_MCoutput_fluorescence.mat'],'I_fluorescence','MCinput_fluorescence','P');
fprintf('./Data/%s_MCoutput_fluorescence.mat saved\n',name);
clear I_fluorescence MCinput_fluorescence

%% Call lookmcxyz
lookmcxyz(name);

return
