function MCoutput = runMonteCarlo(name,varargin)
%
%   Prepares the illumination beam and runs the Monte Carlo simulation.
%   After finishing, calls lookMCmatlab for the display of the result.
%
%   Define the time requested for simulating photons.
%   Define the behaviour of photons that stray outside the cuboid.
%   Define the beam parameters.
%   Depending on the chosen beam type, define the focus position, the beam
%   waist of the focus, the divergence angle and/or the direction of the
%   beam's center axis.
%
%   Input
%       name
%           the basename of the file saved by makeGeometry.m
%       varargin
%           if 'silent' is specified as an additional argument, disables
%           overwrite prompt, command window text and progress indication
%           if 'dontMaxCPU' is specified as an additional argument, ensures
%           that on Windows, MCmatlab will leave one processor unused. (On
%           Mac, multithreading is not yet implemented anyway)
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

%% Check if silent mode or dontMaxCPU was specified
silentMode = double(any(strcmp(varargin,'silent')));
dontMaxCPU = double(any(strcmpi(varargin,'dontMaxCPU')));

%% Check for preexisting files
if(~silentMode && ~deleteDataFiles(name)); return; end

%% Load data from makeGeometry.m
load(['./Data/' name '.mat']);

%% Define parameters (user-specified)
% Simulation duration
simulationTime = .1;      % [min] time duration of the simulation

% Beam type
% 0: Pencil beam
% 1: Isotropically emitting point source
% 2: Infinite plane wave
% 3: Gaussian focus, Gaussian far field beam
% 4: Gaussian focus, top-hat far field beam
% 5: Top-hat focus, Gaussian far field beam
% 6: Top-hat focus, top-hat far field beam
% 7: Laguerre-Gaussian LG01 beam
beamtypeFlag = 1;

% Position of focus, only used for beamtypeFlag ~=2 (if beamtypeFlag == 1 this is the source position)
xFocus = 0;                % [cm] x position of focus
yFocus = 0;                % [cm] y position of focus
zFocus = 0.04;          % [cm] z position of focus

% Direction of beam center axis, only used if beamtypeflag ~= 1:
% Given in terms of the spherical coordinates theta and phi measured in radians, using the ISO
% convention illustrated at https://en.wikipedia.org/wiki/Spherical_coordinate_system
% Keep in mind that the z-axis in the volumetric plots is shown pointing down, so you want to satisfy 0<=theta<pi/2.
% Examples: theta = 0, phi = 0 means a beam going straight down (positive z direction)
%           theta = pi/4, phi = 0 means a beam going halfway between the positive x and positive z directions.
%           theta = pi/4, phi = -pi/2 means a beam going halfway between the negative y and positive z directions
thetaBeam = 0; % [rad]
phiBeam   = 0; % [rad]

% Focus properties and divergence angles, only used if beamtypeflag > 2
waist = 0.03;                  % [cm] focus waist 1/e^2 radius
% divergence = G.wavelength*1e-9/(pi*waist*1e-2); % [rad] Diffraction limited divergence angle for Gaussian beam
divergence = 15/180*pi;         % [rad] divergence 1/e^2 half-angle of beam

%% Optional light collector properties (user-specified)
% A "light collector" in this context can be either an objective lens or a fiber tip
useLightCollector = 0; % Set to 1 for true, 0 for false

% Position of either the center of the objective lens focal plane or the fiber tip
xFPC_LC = 0;% [cm]
yFPC_LC = 0;% [cm]
zFPC_LC = G.nz*G.dz/2;% [cm] negative values correspond to a location above the volume

% Direction that the light collector is facing, defined in the same way as the beam direction using
% ISO spherical coordinates.
% For above the surface (negative z), you'll want to satisfy 0<=theta<=pi/2.
% For below the surface (positive z), such as measuring in transmission, you'll want pi/2<=theta<=pi
% If theta = 0 or pi, phi serves only to rotate the view. If you want the light collector X axis 
% to coincide with the cuboid x axis, use phi = pi/2.
theta_LC = atan(1/sqrt(2)); % [rad]
phi_LC   = -3*pi/4; % [rad]

% Focal length of the objective lens (if light collector is a fiber, set this to Inf).
f_LC = 1; % [cm]

% Diameter of the light collector aperture. For an ideal thin lens, this is 2*f*tan(asin(lensNA)).
diam_LC = 2; % [cm]

% Field Size of the imaging system (diameter of area in object plane that gets imaged). Only used for finite f_LC.
FieldSize_LC = 2; % [cm]

% Fiber NA. Only used for infinite f_LC.
NA_LC = 0.22; % [-]

% Resolution of light collector in pixels
resX_LC = 400;
resY_LC = 400;

%% Check to ensure that the light collector is not inside the cuboid
if useLightCollector
    if isfinite(f_LC)
        xLCC = xFPC_LC - f_LC*sin(theta_LC)*cos(phi_LC); % x position of Light Collector Center
        yLCC = yFPC_LC - f_LC*sin(theta_LC)*sin(phi_LC); % y position
        zLCC = zFPC_LC - f_LC*cos(theta_LC);             % z position
    else
        xLCC = xFPC_LC;
        yLCC = yFPC_LC;
        zLCC = zFPC_LC;
    end

    if (abs(xLCC)               < G.nx*G.dx/2 && ...
        abs(yLCC)               < G.ny*G.dy/2 && ...
        abs(zLCC - G.nz*G.dz/2) < G.nz*G.dz/2)
        error('Error: Light collector center (%.4f,%.4f,%.4f) is inside cuboid',xLCC,yLCC,zLCC);
    end    
end

%% Call Monte Carlo C script (MEX file) to get fluence rate (intensity) distribution
G.M = G.M - 1; % The medium matrix has to be converted from MATLAB's 1-based indexing to C's 0-based indexing
MCinput = struct('silentMode',silentMode,'dontMaxCPU',dontMaxCPU,'G',G,'simulationTime',simulationTime,...
    'beamtypeFlag',beamtypeFlag,'xFocus',xFocus,'yFocus',yFocus,'zFocus',zFocus,'thetaBeam',thetaBeam,...
    'phiBeam',phiBeam,'waist',waist,'divergence',divergence,'useLightCollector',useLightCollector,...
    'xFPC_LC',xFPC_LC,'yFPC_LC',yFPC_LC,'zFPC_LC',zFPC_LC,'theta_LC',theta_LC,'phi_LC',phi_LC,'f_LC',f_LC,...
    'diam_LC',diam_LC,'FieldSize_LC',FieldSize_LC,'NA_LC',NA_LC,'resX_LC',resX_LC,'resY_LC',resY_LC);
MCoutput = MCmatlab(MCinput);

%% Save output and clear memory
save(['./Data/' name '_MCoutput.mat'],'MCoutput','MCinput');
if(~silentMode) fprintf('./Data/%s_MCoutput.mat saved\n',name); end
clear MCoutput

%% Make plots
if(~silentMode) lookMCmatlab(name); end

end
