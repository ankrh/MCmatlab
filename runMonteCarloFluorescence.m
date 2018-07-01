function runMonteCarloFluorescence(name)
%   Script for simulating distribution and magnitude of fluorescence
%   based on the output of runMonteCarlo.m
%
%   Prepares and runs the Monte Carlo simulation.
%   After finishing, calls lookMCmatlab for the display of the result.
%
%   Define the time requested for simulating photons.
%   Define the power of the excitation beam (used to take saturation into
%   account).
%   Define the behaviour of photons that stray outside the tissue cuboid:
%       0 = no boundaries: photons wander freely also outside the tissue
%       cuboid and get killed only if they wander too far (6 times the cuboid
%       size).
%       1 = escape at boundaries: photons that stray outside the tissue
%       cuboid get killed immediately.
%       2 = escape at surface only: photons that hit the top surface get
%       killed immediately, photons hitting other surfaces can wander up to
%       6 times the cuboid size.
%
%   Input
%       name
%           the basename of the file as specified in makeTissue.m
%
%   Output
%       ./Data/[name]_MCoutput_fluorescence.mat
%           file containing the 3D fluorescence fluence rate distribution
%
%   Requires
%       deleteDataFiles.m
%       MCmatlab.mex (architecture specific)
%       lookMCmatlab.m
%

%% Updates
% Created 2018-04-19 by Anders K. Hansen, DTU Fotonik

%% Check for preexisting files
if(~deleteDataFiles(name)); return; end

%% Load data from makeTissue.m and runMonteCarlo.m
load(['./Data/' name '.mat']);
load(['./Data/' name '_MCoutput.mat']);

MCinput_fluorescence = MCinput;

%% Define parameters (user-specified)
% Simulation duration
MCinput_fluorescence.simulationTime = 0.1;      % [min] time duration of the simulation

% Incident excitation power, for taking fluorescence saturations into account
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
MCinput_fluorescence.sourceDistribution = Y_vec(T).*mua_vec(T)*P.*MCoutput.F./(1 + P*MCoutput.F./sat_vec(T)); % [W/cm^3]

%% Call Monte Carlo C script to get fluorescence distribution
MCoutput_fluorescence = MCmatlab(MCinput_fluorescence); % MCoutput_fluorescence.F is an absolute fluence rate (intensity) quantity, unlike the non-fluorescence MCoutput.F which are actually fluence rates normalized to the incident power

%% Save output and clear memory
save(['./Data/' name '_MCoutput_fluorescence.mat'],'MCoutput_fluorescence','MCinput_fluorescence','P');
fprintf('./Data/%s_MCoutput_fluorescence.mat saved\n',name);
clear MCoutput_fluorescence MCinput_fluorescence

%% Make plots
lookMCmatlab(name);

end
