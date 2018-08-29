%% Geometry definition
clear Ginput
Ginput.silentMode = false;
Ginput.assumeMatchedInterfaces = true;
Ginput.boundaryType = 1;

Ginput.wavelength  = 532;		% [nm] Wavelength of the Monte Carlo simulation

Ginput.nx = 100;				% number of bins in the x direction
Ginput.ny = 100;				% number of bins in the y direction
Ginput.nz = 100;				% number of bins in the z direction
Ginput.Lx = .1;				% [cm] x size of simulation area
Ginput.Ly = .1;				% [cm] y size of simulation area
Ginput.Lz = .1;				% [cm] z size of simulation area

Ginput.GeomFunc = @GeometryDefinition_BloodVessel; % Specify which function (defined at the end of this m file) to use for defining the distribution of media in the cuboid

% Execution, do not modify the next line:
G = defineGeometry(Ginput);

%% Monte Carlo simulation
clear MCinput
MCinput.silentMode = false;
MCinput.useAllCPUs = true;
MCinput.simulationTime = .1;      % [min] time duration of the simulation

MCinput.Beam.beamType = 2;
MCinput.Beam.xFocus = 0;                % [cm] x position of focus
MCinput.Beam.yFocus = 0;                % [cm] y position of focus
MCinput.Beam.zFocus = Ginput.Lz/2;      % [cm] z position of focus
MCinput.Beam.theta = 0; % [rad]
MCinput.Beam.phi   = 0; % [rad]
MCinput.Beam.waist = 0.005;                  % [cm] focus waist 1/e^2 radius
MCinput.Beam.divergence = 5/180*pi;         % [rad] divergence 1/e^2 half-angle of beam (for a diffraction limited Gaussian beam, this is G.wavelength*1e-9/(pi*MCinput.Beam.waist*1e-2))

% Execution, do not modify the next two lines:
MCinput.G = G;
MCoutput = runMonteCarlo(MCinput);

%% Explanations
% silentMode = true disables overwrite prompt,
% command window text, progress indication and plot generation

% If assumeMatchedInterfaces = true, all refractive indices
% are assumed to be 1 and there is no Fresnel reflection or refraction.
% Otherwise, refractive indices from getMediaProperties are used. Note that
% non-matched interfaces must be normal to the z axis, so each xy-slice
% must have a constant refractive index. 

% Boundary type
% 0: No boundaries. Photons are allowed to leave the cuboid and are still
%    tracked outside, including absorption and scattering events. They get
%    terminated only if they wander too far (6 times the cuboid size).
% 1: Cuboid boundaries. All 6 cuboid surfaces are considered photon boundaries.
% 2: Top boundary only. Only the top surface (z = 0) is a photon boundary.
% Regardless of the boundary type, photons that wander 6 times the cuboid
% size will be terminated. When a photon hits a photon boundary at a position
% where the refractive index is 1, it escapes and may contribute to the
% signal of the light collector depending on its trajectory. Otherwise, the
% photon is just terminated, meaning that it cannot contribute to the light
% collector.

% If useAllCPUs = true, MCmatlab will use all available processors on Windows. Otherwise,
% one will be left unused. Useful for doing other work on the PC
% while simulations are running.

% Beam type
% 0: Pencil beam
% 1: Isotropically emitting point source
% 2: Infinite plane wave
% 3: Gaussian focus, Gaussian far field beam
% 4: Gaussian focus, top-hat far field beam
% 5: Top-hat focus, Gaussian far field beam
% 6: Top-hat focus, top-hat far field beam
% 7: Laguerre-Gaussian LG01 beam

% xFocus, yFocus, zFocus: Position of focus in the absence of any refraction, only used for beamType ~=2 (if beamType == 1 this is the source position)

% theta and phi define the direction of beam center axis, only used if beamtypeflag ~= 1:
% Given in terms of the spherical coordinates theta and phi measured in radians, using the ISO
% convention illustrated at https://en.wikipedia.org/wiki/Spherical_coordinate_system
% Keep in mind that the z-axis in the volumetric plots is shown pointing down, so you want to satisfy 0<=theta<pi/2.
% Examples: theta = 0, phi = 0 means a beam going straight down (positive z direction)
%           theta = pi/4, phi = 0 means a beam going halfway between the positive x and positive z directions.
%           theta = pi/4, phi = -pi/2 means a beam going halfway between the negative y and positive z directions



function M = GeometryDefinition_BloodVessel(X,Y,Z,parameters)
% Blood vessel example:
zsurf = 0.01;
epd_thick = 0.006;
vesselradius  = 0.0100;
vesseldepth = 0.04;
M = 2*ones(size(X)); % fill background with water (gel)
M(Z > zsurf) = 4; % epidermis
M(Z > zsurf + epd_thick) = 5; % dermis
M(X.^2 + (Z - (zsurf + vesseldepth)).^2 < vesselradius^2) = 6; % blood
end
