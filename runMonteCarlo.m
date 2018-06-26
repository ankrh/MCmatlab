function runMonteCarlo(name)
%
%   Prepares the illumination beam and runs the Monte Carlo simulation.
%   After finishing, calls lookmcxyz for the display of the result.
%
%   Define the time requested for simulating photons.
%   Define the behaviour of photons that stray outside the tissue cuboid:
%       0 = no boundaries: photons wander freely also outside the tissue
%       cuboid and get killed only if they wander too far (6 times the cuboid
%       size).
%       1 = escape at boundaries: photons that stray outside the tissue
%       cuboid get killed immediately.
%       2 = escape at surface only: photons that hit the top surface get
%       killed immediately, photons hitting other surfaces can wander up to
%       6 times the cuboid size.
%   Define the beam parameters. The following beam types can be requested:
%       0 = top-hat focus, top-hat far field beam
%       1 = Gaussian focus, Gaussian far field beam
%       2 = isotropically emitting point
%       3 = infinite plane wave
%       4 = pencil beam
%       5 = top-hat focus, Gaussian far field beam
%       6 = Gaussian focus, top-hat far field beam
%       7 = Laguerre-Gaussian LG01 focus, LG01 far field beam
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
%       mcxyz.mex (architecture specific)
%       lookmcxyz.m
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
MCinput.simulationTime = .1;      % [min] time duration of the simulation

% Boundary type
% 0 = no boundaries
% 1 = escape at boundaries
% 2 = escape at surface only. No x, y, bottom z boundaries
MCinput.boundaryFlag = 1;

% Beam type
% 0 = top-hat focus, top-hat far field beam
% 1 = Gaussian focus, Gaussian far field beam
% 2 = isotropically emitting point
% 3 = infinite plane wave
% 4 = pencil beam
% 5 = top-hat focus, Gaussian far field beam
% 6 = Gaussian focus, top-hat far field beam
% 7 = Laguerre-Gaussian LG01 focus, LG01 far field beam
MCinput.beamtypeFlag = 3;

% Position of focus, only used for beamtypeFlag ~=3 (if beamtypeFlag == 2 this is the source position)
MCinput.xFocus = 0;                % [cm] x position of focus
MCinput.yFocus = 0;                % [cm] y position of focus
MCinput.zFocus = nz*dz/2;          % [cm] z position of focus

% Direction of beam center axis, only used if beamtypeflag ~= 2
% Given in terms of spherical coordinates theta and phi measured in radians, using the ISO
% convention illustrated at https://en.wikipedia.org/wiki/Spherical_coordinate_system
% Keep in mind that th z-axis in the volumetric plots is shown pointing down.
% For basically all cases, you want to satisfy 0<=theta<pi/2.
% Examples: theta = 0, phi = 0 means a beam going straight down (positive z direction)
%           theta = pi/4, phi = 0 means a beam going halfway between the positive x and positive z directions.
%           theta = pi/4, phi = -pi/2 means a beam going halfway between the negative y and positive z directions
MCinput.thetaBeam = 0; % [rad]
MCinput.phiBeam   = 0; % [rad]

% Focus properties and divergence angles, only used if beamtypeflag == 0, 1, 5, 6 or 7
MCinput.waist = 0.015;             % [cm] focus waist 1/e^2 radius
% MCinput.divergence = wavelength*1e-9/(pi*MCinput.waist*1e-2); % [rad] Diffraction limited divergence angle for Gaussian beam
MCinput.divergence = 15/180*pi;         % [rad] divergence 1/e^2 half-angle of beam

%% Optional detector properties
% Position of either the center of the objective lens focal plane or the fiber tip
MCinput.xFPCDet = 0; % [cm]
MCinput.yFPCDet = 0; % [cm]
MCinput.zFPCDet = 0; % [cm] negative values correspond to a location above the volume

% Direction that the detector is facing, defined in the same way as the beam direction using
% ISO spherical coordinates
% For detector above the surface, you usually want to satisfy 0<=theta<=pi/2.
% For a detector below the surface, such as measuring in transmission, you usually want pi/2<=theta<=pi
% If theta = 0 or pi, phi serves only to rotate the view of the tissue. If you want the detector X axis 
% to coincide with the tissue cuboid x axis, use phi = pi/2 in that case.
MCinput.thetaDet = 0; % [rad]
MCinput.phiDet   = pi/2; % [rad]

% Focal length of the detector (if it is not a lens, set this to inf). Must be positive.
MCinput.fDet = inf; % [cm]

% Diameter of the detector aperture. For an ideal thin lens, this is 2*tan(arcsin(lensNA/f)).
MCinput.diamDet = 0.5; % [cm]

% For an objective lens: Field Size of the imaging system (diameter of area in object plane that gets imaged). For a fiber tip: The fiber's NA.
MCinput.FSorNADet = 2; % [cm or dimensionless]

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
[F,LCP,ImP,LCFF] = mcxyz(MCinput);

%% Save output and clear memory
save(['./Data/' name '_MCoutput.mat'],'F','MCinput');
fprintf('./Data/%s_MCoutput.mat saved\n',name);
% clear F MCinput

%% Make plots
close all
lookmcxyz(name);

figure;imagesc([-MCinput.diamDet MCinput.diamDet]/2,[-MCinput.diamDet MCinput.diamDet]/2,LCP.');title('Light Collector Plane');axis xy;axis equal;axis tight;xlabel('X [cm]');ylabel('Y [cm]');
figure;imagesc([-MCinput.FSorNADet MCinput.FSorNADet]/2,[-MCinput.FSorNADet MCinput.FSorNADet]/2,ImP.');title('Image Plane');axis xy;axis equal;axis tight;xlabel('X [cm]');ylabel('Y [cm]');
figure;imagesc([0 pi/2],[-pi pi],LCFF.');title('Light Collector Far Field');axis xy;axis equal;axis tight;xlabel('theta [rad]');ylabel('phi [rad]');

end
