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

%% Optional light collector properties (user-specified)
MCinput_fluorescence.useLightCollector = 1; % Set to 1 for true, 0 for false

% Position of either the center of the objective lens focal plane or the fiber tip
MCinput_fluorescence.xFPC_LC = 0; % [cm]
MCinput_fluorescence.yFPC_LC = 0; % [cm]
MCinput_fluorescence.zFPC_LC = nz*dz/2; % [cm] negative values correspond to a location above the volume

% Direction that the light collector is facing, defined in the same way as the beam direction using
% ISO spherical coordinates.
% For above the surface (negative z), you'll want to satisfy 0<=theta<=pi/2.
% For below the surface (positive z), such as measuring in transmission, you'll want pi/2<=theta<=pi
% If theta = 0 or pi, phi serves only to rotate the view of the tissue. If you want the light collector X axis 
% to coincide with the tissue cuboid x axis, use phi = pi/2.
MCinput_fluorescence.theta_LC = 0; % [rad]
MCinput_fluorescence.phi_LC   = pi/2; % [rad]

% Focal length of the objective lens (if light collector is a fiber, set this to inf).
MCinput_fluorescence.f_LC = 1; % [cm]

% Diameter of the light collector aperture. For an ideal thin lens, this is 2*f*tan(asin(lensNA)).
MCinput_fluorescence.diam_LC = 2; % [cm]

% Field Size of the imaging system (diameter of area in object plane that gets imaged). Only used for finite f_LC.
MCinput_fluorescence.FieldSize_LC = 2; % [cm]

% Fiber NA. Only used for infinite f_LC.
MCinput_fluorescence.NA_LC = 0.22; % [-]

% Resolution of light collector in pixels
MCinput_fluorescence.resX_LC = 200;
MCinput_fluorescence.resY_LC = 200;

%% Set remaining parameters
MCinput_fluorescence.tissueList = tissueList_fluorescence;

%% Calculate 3D fluorescence source distribution, including saturation
mua_vec = [tissueList.mua]; % The tissues' excitation absorption coefficients
Y_vec = [tissueList.Y]; % The tissues' fluorescence power yields
sat_vec = [tissueList.sat]; % The tissues' fluorescence saturation fluence rates (intensity)
MCinput_fluorescence.sourceDistribution = Y_vec(T).*mua_vec(T)*P.*MCoutput.F./(1 + P*MCoutput.F./sat_vec(T)); % [W/cm^3]

%% Check to ensure that the light collector is not inside the cuboid
if MCinput_fluorescence.useLightCollector
    if isfinite(MCinput_fluorescence.f_LC)
        xLCC = MCinput_fluorescence.xFPC_LC - MCinput_fluorescence.f_LC*sin(MCinput_fluorescence.theta_LC)*cos(MCinput_fluorescence.phi_LC); % x position of Light Collector Center
        yLCC = MCinput_fluorescence.yFPC_LC - MCinput_fluorescence.f_LC*sin(MCinput_fluorescence.theta_LC)*sin(MCinput_fluorescence.phi_LC); % y position
        zLCC = MCinput_fluorescence.zFPC_LC - MCinput_fluorescence.f_LC*cos(MCinput_fluorescence.theta_LC);
    else
        xLCC = MCinput_fluorescence.xFPC_LC;
        yLCC = MCinput_fluorescence.yFPC_LC;
        zLCC = MCinput_fluorescence.zFPC_LC;
    end

    if (abs(xLCC)           < nx*dx/2 && ...
        abs(yLCC)           < ny*dy/2 && ...
        abs(zLCC - nz*dz/2) < nz*dz/2)
        error('Error: Light collector center (%.4f,%.4f,%.4f) is inside cuboid',xLCC,yLCC,zLCC);
    end    
end

%% Call Monte Carlo C script to get fluorescence distribution
MCoutput_fluorescence = MCmatlab(MCinput_fluorescence); % MCoutput_fluorescence.F is an absolute fluence rate (intensity) quantity, unlike the non-fluorescence MCoutput.F which are actually fluence rates normalized to the incident power

%% Save output and clear memory
save(['./Data/' name '_MCoutput_fluorescence.mat'],'MCoutput_fluorescence','MCinput_fluorescence','P');
fprintf('./Data/%s_MCoutput_fluorescence.mat saved\n',name);
clear MCoutput_fluorescence MCinput_fluorescence

%% Make plots
lookMCmatlab(name);

end
