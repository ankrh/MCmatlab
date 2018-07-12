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
MCinput.boundaryFlag = 1;

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
MCinput.zFocus = 0.028;          % [cm] z position of focus

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
MCinput.waist = 0.000003;                  % [cm] focus waist 1/e^2 radius
% MCinput.divergence = wavelength*1e-9/(pi*MCinput.waist*1e-2); % [rad] Diffraction limited divergence angle for Gaussian beam
MCinput.divergence = 25/180*pi;         % [rad] divergence 1/e^2 half-angle of beam

%% Optional light collector properties (user-specified)
% A "light collector" in this context can be either an objective lens or a fiber tip
MCinput.useLightCollector = 0; % Set to 1 for true, 0 for false

% Position of either the center of the objective lens focal plane or the fiber tip
MCinput.xFPC_LC = 0;% [cm]
MCinput.yFPC_LC = 0;% [cm]
MCinput.zFPC_LC = nz*dz/2;% [cm] negative values correspond to a location above the volume

% Direction that the light collector is facing, defined in the same way as the beam direction using
% ISO spherical coordinates.
% For above the surface (negative z), you'll want to satisfy 0<=theta<=pi/2.
% For below the surface (positive z), such as measuring in transmission, you'll want pi/2<=theta<=pi
% If theta = 0 or pi, phi serves only to rotate the view of the tissue. If you want the light collector X axis 
% to coincide with the tissue cuboid x axis, use phi = pi/2.
MCinput.theta_LC = atan(1/sqrt(2)); % [rad]
MCinput.phi_LC   = -3*pi/4; % [rad]

% Focal length of the objective lens (if light collector is a fiber, set this to inf).
MCinput.f_LC = 1; % [cm]

% Diameter of the light collector aperture. For an ideal thin lens, this is 2*f*tan(asin(lensNA)).
MCinput.diam_LC = 2; % [cm]

% Field Size of the imaging system (diameter of area in object plane that gets imaged). Only used for finite f_LC.
MCinput.FieldSize_LC = 2; % [cm]

% Fiber NA. Only used for infinite f_LC.
MCinput.NA_LC = 0.22; % [-]

% Resolution of light collector in pixels
MCinput.resX_LC = 400;
MCinput.resY_LC = 400;

%% Determine remaining parameters
% Voxel sizes
MCinput.dx = x(2)-x(1);            % [cm] voxel size in x direction
MCinput.dy = y(2)-y(1);            % [cm] voxel size in y direction
MCinput.dz = z(2)-z(1);            % [cm] voxel size in z direction

% Tissue definition
MCinput.tissueList = tissueList;
MCinput.T = T-1; % The tissue matrix has to be converted from MATLAB's 1-based indexing to C's 0-based indexing
clear T

% Array of refractive indices along z axis
MCinput.RI = RI;

%% Check to ensure that the light collector is not inside the cuboid
if MCinput.useLightCollector
    if isfinite(MCinput.f_LC)
        xLCC = MCinput.xFPC_LC - MCinput.f_LC*sin(MCinput.theta_LC)*cos(MCinput.phi_LC); % x position of Light Collector Center
        yLCC = MCinput.yFPC_LC - MCinput.f_LC*sin(MCinput.theta_LC)*sin(MCinput.phi_LC); % y position
        zLCC = MCinput.zFPC_LC - MCinput.f_LC*cos(MCinput.theta_LC);
    else
        xLCC = MCinput.xFPC_LC;
        yLCC = MCinput.yFPC_LC;
        zLCC = MCinput.zFPC_LC;
    end

    if (abs(xLCC)           < nx*dx/2 && ...
        abs(yLCC)           < ny*dy/2 && ...
        abs(zLCC - nz*dz/2) < nz*dz/2)
        error('Error: Light collector center (%.4f,%.4f,%.4f) is inside cuboid',xLCC,yLCC,zLCC);
    end    
end

%% Call Monte Carlo C script (mex file) to get fluence rate (intensity) distribution
MCoutput = MCmatlab(MCinput);

%% Save output and clear memory
save(['./Data/' name '_MCoutput.mat'],'MCoutput','MCinput');
fprintf('./Data/%s_MCoutput.mat saved\n',name);
clear MCoutput MCinput

%% Make plots
lookMCmatlab(name);

end
