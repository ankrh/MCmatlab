function MCoutput_fluorescence = runMonteCarloFluorescence(name,varargin)
%
%   Script for simulating distribution and magnitude of fluorescence
%   based on the output of runMonteCarlo.m
%
%   Prepares and runs the Monte Carlo simulation.
%   After finishing, calls lookMCmatlab for the display of the result.
%
%   Define the time requested for simulating photons.
%   Define the power of the excitation beam (used to take saturation into
%   account).
%   Define the behaviour of photons that stray outside the cuboid.
%
%   Input
%       name
%           the basename of the files saved by makeGeometry.m and runMonteCarlo.m
%       varargin
%           if 'silent' is specified as an additional argument, disables
%           overwrite prompt, command window text and progress indication
%           if 'dontMaxCPU' is specified as an additional argument, ensures
%           that on Windows, MCmatlab will leave one processor unused. (On
%           Mac, multithreading is not yet implemented anyway)
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

%% Load data from makeGeometry.m and runMonteCarlo.m
load(['./Data/' name '.mat']);
load(['./Data/' name '_MCoutput.mat'],'MCoutput');

%% Check if silent mode or dontMaxCPU was specified
silentMode = double(any(strcmp(varargin,'silent')));
dontMaxCPU = double(any(strcmpi(varargin,'dontMaxCPU')));

%% Check for preexisting files
if(~silentMode && ~deleteDataFiles(name)); return; end

%% Define parameters (user-specified)
% Simulation duration
simulationTime = 0.1;      % [min] time duration of the simulation

% Incident excitation power, for taking fluorescence saturations into account
P = 2; % [W]

%% Optional light collector properties (user-specified)
% A "light collector" in this context can be either an objective lens or a fiber tip
useLightCollector = 1; % Set to 1 for true, 0 for false

% Position of either the center of the objective lens focal plane or the fiber tip
xFPC_LC = 0; % [cm]
yFPC_LC = 0; % [cm]
zFPC_LC = G.nz*G.dz/2; % [cm] negative values correspond to a location above the volume

% Direction that the light collector is facing, defined in the same way as the beam direction using
% ISO spherical coordinates.
% For above the surface (negative z), you'll want to satisfy 0<=theta<=pi/2.
% For below the surface (positive z), such as measuring in transmission, you'll want pi/2<=theta<=pi
% If theta = 0 or pi, phi serves only to rotate the view. If you want the light collector X axis 
% to coincide with the cuboid x axis, use phi = pi/2.
theta_LC = 0; % [rad]
phi_LC   = pi/2; % [rad]

% Focal length of the objective lens (if light collector is a fiber, set this to Inf).
f_LC = 1; % [cm]

% Diameter of the light collector aperture. For an ideal thin lens, this is 2*f*tan(asin(lensNA)).
diam_LC = 2; % [cm]

% Field Size of the imaging system (diameter of area in object plane that gets imaged). Only used for finite f_LC.
FieldSize_LC = 2; % [cm]

% Fiber NA. Only used for infinite f_LC.
NA_LC = 0.22; % [-]

% Resolution of light collector in pixels
resX_LC = 200;
resY_LC = 200;

%% Set remaining parameters
mediaProperties = mediaProperties_fluorescence;

%% Calculate 3D fluorescence source distribution, including saturation
mua_vec = [mediaProperties_fluorescence.mua]; % The media's excitation absorption coefficients
Y_vec = [mediaProperties_fluorescence.Y]; % The media's fluorescence power yields
sat_vec = [mediaProperties_fluorescence.sat]; % The media's fluorescence saturation fluence rates (intensity)
sourceDistribution = Y_vec(T).*mua_vec(T)*P.*MCoutput.F./(1 + P*MCoutput.F./sat_vec(T)); % [W/cm^3]

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

%% Call Monte Carlo C script (MEX file) to get fluorescence fluence rate (intensity) distribution
G.M = G.M - 1; % The medium matrix has to be converted from MATLAB's 1-based indexing to C's 0-based indexing
MCinput = struct('silentMode',silentMode,'dontMaxCPU',dontMaxCPU,'G',G,'simulationTime',simulationTime,...
    'sourceDistribution',sourceDistribution,'useLightCollector',useLightCollector,...
    'xFPC_LC',xFPC_LC,'yFPC_LC',yFPC_LC,'zFPC_LC',zFPC_LC,'theta_LC',theta_LC,'phi_LC',phi_LC,'f_LC',f_LC,...
    'diam_LC',diam_LC,'FieldSize_LC',FieldSize_LC,'NA_LC',NA_LC,'resX_LC',resX_LC,'resY_LC',resY_LC);
MCoutput_fluorescence = MCmatlab(MCinput_fluorescence); % MCoutput_fluorescence.F is an absolute fluence rate (intensity) quantity, unlike the non-fluorescence MCoutput.F which are actually fluence rates normalized to the incident power

%% Save output and clear memory
save(['./Data/' name '_MCoutput_fluorescence.mat'],'MCoutput_fluorescence','MCinput_fluorescence','P');
if(~silentMode) fprintf('./Data/%s_MCoutput_fluorescence.mat saved\n',name); end
clear MCoutput_fluorescence

%% Make plots
if(~silentMode) lookMCmatlab(name); end

end
