function runMonteCarlo(name)
%
%   Prepares the illumination beam and runs the Monte Carlo simulation.
%   After finishing, calls lookMCmatlab for the display of the result.
%
%   Define the time requested for simulating photons.
%   Define the behaviour of photons that stray outside the tissue cuboid.
%   Define the beam parameters.
%   Depending on the chosen beam type, define the focus position, the beam
%   waist of the focus, the divergence angle and/or the direction of the
%   beam's center axis.
%
%   Input
%       name
%           the basename of the file as specified in makeTissue.m
%
%   Output
%       ./Data/[name]_MCoutput.mat
%           file containing the 3D fluence rate distribution
%
%   Requires
%       deleteDataFiles.m
%       MCmatlab.mex (architecture specific)
%       lookMCmatlab.m
%

%% Updates
%   2014-01: Mathias Christensen & Rasmus L. Pedersen, DTU Fotonik
%   2017-06: Anders K. Hansen & Dominik Marti, DTU Fotonik
%   2018-04: Anders K. Hansen

%% Check for preexisting files
if(~deleteDataFiles(name)); return; end

%% Load data from makeTissue.m
load(['./Data/' name '.mat']);

%% Define parameters (user-specified)
% Simulation duration
MCinput.simulationTime = 0.1;      % [min] time duration of the simulation

% Boundary type
% 0: No boundaries. Photons wander freely also outside the tissue cuboid and
%    get killed only if they wander too far (6 times the cuboid size).
% 1: Escape at boundaries. Photons that stray outside the tissue cuboid get
%    killed immediately.
% 2: Escape at surface only. Photons that hit the top surface get killed
%    immediately, photons hitting other surfaces can wander up to 6 times
%    the cuboid size.
MCinput.boundaryFlag = 2;

% Beam type
% 0: Pencil beam
% 1: Isotropically emitting point source
% 2: Infinite plane wave
% 3: Gaussian focus, Gaussian far field beam
% 4: Gaussian focus, top-hat far field beam
% 5: Top-hat focus, Gaussian far field beam
% 6: Top-hat focus, top-hat far field beam
% 7: Laguerre-Gaussian LG01 beam
MCinput.beamtypeFlag = 6;

% Position of focus, only used for beamtypeFlag ~=2 (if beamtypeFlag == 1 this is the source position)
MCinput.xFocus = 0;                % [cm] x position of focus
MCinput.yFocus = 0;                % [cm] y position of focus
MCinput.zFocus = 0;                % [cm] z position of focus

% Direction of beam center axis, only used if beamtypeflag ~= 1:
% Given in terms of the spherical coordinates theta and phi measured in radians, using the ISO
% convention illustrated at https://en.wikipedia.org/wiki/Spherical_coordinate_system
% Keep in mind that the z-axis in the volumetric plots is shown pointing down, so you want to satisfy 0<=theta<pi/2.
% Examples: theta = 0, phi = 0 means a beam going straight down (positive z direction)
%           theta = pi/4, phi = 0 means a beam going halfway between the positive x and positive z directions.
%           theta = pi/4, phi = -pi/2 means a beam going halfway between the negative y and positive z directions
MCinput.thetaBeam = 0; % [rad]
MCinput.phiBeam   = 0; % [rad]

% Focus properties and divergence angles, only used if beamtypeflag > 2
MCinput.waist = 0.03;             % [cm] focus waist 1/e^2 radius
% MCinput.divergence = wavelength*1e-9/(pi*MCinput.waist*1e-2); % [rad] Diffraction limited divergence angle for Gaussian beam
MCinput.divergence = 0/180*pi;         % [rad] divergence 1/e^2 half-angle of beam

%% Determine remaining parameters
% Voxel sizes
MCinput.dx = x(2)-x(1);            % [cm] voxel size in x direction
MCinput.dy = y(2)-y(1);            % [cm] voxel size in y direction
MCinput.dz = z(2)-z(1);            % [cm] voxel size in z direction

% Tissue definition
MCinput.tissueList = tissueList;
MCinput.T = T-1; % The tissue matrix has to be converted from MATLAB's 1-based indexing to C's 0-based indexing
clear T

%% Call Monte Carlo C script (mex file) to get fluence rate (intensity) distribution
F = MCmatlab(MCinput);

%% Save output and clear memory
save(['./Data/' name '_MCoutput.mat'],'F','MCinput');
fprintf('./Data/%s_MCoutput.mat saved\n',name);
clear F MCinput

%% Make plots
lookMCmatlab(name);

end
