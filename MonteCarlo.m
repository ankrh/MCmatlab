function MonteCarlo(name)
%% Load data from makeTissue.m
load(['.\Data\' name '.mat']);

%% Define parameters (user-specified)
% Simulation duration
MCinput.simulationTime = 0.5;      % [min] time duration of the simulation

% Beam type
MCinput.beamtypeFlag = 3;          % beam type: 0 = top-hat focus, top-hat far field beam, 1 = Gaussian focus, Gaussian far field beam, 2 = isotropically emitting point, 3 = infinite plane wave, 4 = pencil beam, 5 = top-hat focus, Gaussian far field beam, 6 = Gaussian focus, top-hat far field beam

% Boundary type
MCinput.boundaryFlag = 1;          % boundary type: 0 = no boundaries, 1 = escape at boundaries, 2 = escape at surface only. No x, y, bottom z boundaries

% Position of focus, only used for beamtypeflag ~=3 (if beamtypeflag == 2 this is the source position)
MCinput.xFocus = 0;                % [cm] x position of focus
MCinput.yFocus = 0;                % [cm] y position of focus
MCinput.zFocus = 0.5;              % [cm] z position of focus

% Direction of beam center axis, only used if beamtypeflag ~= 2:
MCinput.ux0 = 0;                   % trajectory unit vector x composant
MCinput.uy0 = 0;                   % trajectory unit vector y composant
MCinput.uz0 = sqrt(1-MCinput.ux0^2-MCinput.uy0^2); % % trajectory unit vector z composant. Make sure that ux0^2 + uy0^2 + uz0^2 = 1.

% Focus properties and divergence angles, only used if beamtypeflag == 0, 1, 5 or 6
MCinput.waist = 0.025;             % [cm] focus waist 1/e^2 radius
MCinput.divergence = pi/8;         % [rad] divergence 1/e^2 half-angle of beam
% MCinput.divergence  = wavelength*1e-9/(pi*waist*1e-2); % [rad] Diffraction limited divergence angle for Gaussian beam

%% Determine remaining parameters
% Voxel sizes
MCinput.dx = x(2)-x(1);            % [cm] voxel size in x direction
MCinput.dy = y(2)-y(1);            % [cm] voxel size in y direction
MCinput.dz = z(2)-z(1);            % [cm] voxel size in z direction

% Tissue definition
MCinput.tissueList = tissueList;
MCinput.T = T-1; % The tissue matrix has to be converted from MATLAB's 1-based indexing to C's 0-based indexing

%% Call Monte Carlo C script (mex file) to get fluence rate (intensity) distribution
F = mcxyz_mex(MCinput);

%% Calculate normalized volumetric power (NVP)
mua_vec = [tissueList.mua];
NVP = mua_vec(T).*F;

%% Save output
save(['.\Data\' name '_MCoutput.mat'],'F','NVP','MCinput');
fprintf('.\\Data\\%s_MCoutput.mat saved\n',name);

%% Call lookmcxyz
lookmcxyz(name);
