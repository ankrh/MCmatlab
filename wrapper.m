%% Geometry definition
Ginput.silentMode = true;
Ginput.assumeMatchedInterfaces = true;
Ginput.boundaryType = 1;

Ginput.wavelength  = 450;		% [nm] Wavelength of the Monte Carlo simulation
Ginput.wavelength_f = 550;		% [nm] Fluorescence wavelength (set this to NaN for simulations without fluorescence)

Ginput.nx = 100;				% number of bins in the x direction
Ginput.ny = 100;				% number of bins in the y direction
Ginput.nz = 100;				% number of bins in the z direction
Ginput.Lx = .1;				% [cm] x size of simulation area
Ginput.Ly = .1;				% [cm] y size of simulation area
Ginput.Lz = .1;				% [cm] z size of simulation area

Ginput.GeomFunc = @GeometryDefinition_FluorescingCylinder; % Specify which function (defined at the end of this m file) to use for defining the distribution of media in the cuboid
% Ginput.GeomFuncParams = {0.03}; % Cell array containing any additional parameters to pass into the geometry function, such as media depths, inhomogeneity positions, radii etc.

%% Monte Carlo simulation
MCinput.silentMode = true;
MCinput.useAllCPUs = true;
MCinput.simulationTime = .5;      % [min] time duration of the simulation

MCinput.Beam.beamType = 2;
MCinput.Beam.xFocus = 0;                % [cm] x position of focus
MCinput.Beam.yFocus = 0;                % [cm] y position of focus
MCinput.Beam.zFocus = Ginput.Lz/2;      % [cm] z position of focus
MCinput.Beam.theta = 0; % [rad]
MCinput.Beam.phi   = 0; % [rad]
MCinput.Beam.waist = 0.005;                  % [cm] focus waist 1/e^2 radius
MCinput.Beam.divergence = 5/180*pi;         % [rad] divergence 1/e^2 half-angle of beam (for a diffraction limited Gaussian beam, this is G.wavelength*1e-9/(pi*MCinput.Beam.waist*1e-2))

MCinput.LightCollector.xFPC = 0; % [cm] x position of either the center of the objective lens focal plane or the fiber tip
MCinput.LightCollector.yFPC = 0; % [cm]
MCinput.LightCollector.zFPC = 0.03; % [cm]

MCinput.LightCollector.theta = 0; % [rad]
MCinput.LightCollector.phi   = pi/2; % [rad]

MCinput.LightCollector.f = .2; % [cm] Focal length of the objective lens (if light collector is a fiber, set this to Inf).
MCinput.LightCollector.diam = .1; % [cm] Diameter of the light collector aperture. For an ideal thin lens, this is 2*f*tan(asin(lensNA)).
MCinput.LightCollector.FieldSize = .1; % [cm] Field Size of the imaging system (diameter of area in object plane that gets imaged). Only used for finite f.
MCinput.LightCollector.NA = 0.22; % [-] Fiber NA. Only used for infinite f.

MCinput.LightCollector.res = 50; % X and Y resolution of light collector in pixels, only used for finite f

% MCinput.LightCollector.tStart = -1e-13; % [s] Start of the detection time interval
% MCinput.LightCollector.tEnd   = 5e-12; % [s] End of the detection time interval
% MCinput.LightCollector.nTimeBins = 30; % Number of bins between tStart and tEnd

%% Fluorescence Monte Carlo
FMCinput.silentMode = false;
FMCinput.useAllCPUs = true;
FMCinput.simulationTime = .5;      % [min] time duration of the simulation

FMCinput.Beam.P_excitation = 2; % [W]

FMCinput.LightCollector.xFPC = 0; % [cm] x position of either the center of the objective lens focal plane or the fiber tip
FMCinput.LightCollector.yFPC = 0; % [cm]
FMCinput.LightCollector.zFPC = 0.03; % [cm]

FMCinput.LightCollector.theta = 0; % [rad]
FMCinput.LightCollector.phi   = pi/2; % [rad]

FMCinput.LightCollector.f = .2; % [cm] Focal length of the objective lens (if light collector is a fiber, set this to Inf).
FMCinput.LightCollector.diam = .1; % [cm] Diameter of the light collector aperture. For an ideal thin lens, this is 2*f*tan(asin(lensNA)).
FMCinput.LightCollector.FieldSize = .1; % [cm] Field Size of the imaging system (diameter of area in object plane that gets imaged). Only used for finite f.
FMCinput.LightCollector.NA = 0.22; % [-] Fiber NA. Only used for infinite f.

FMCinput.LightCollector.res = 50; % X and Y resolution of light collector in pixels, only used for finite f

%% Execution, do not modify this
MCinput.G = defineGeometry(Ginput);
clear Ginput
FMCinput.MCoutput = runMonteCarlo(MCinput);
FMCinput.G = MCinput.G;
clear MCinput
FMCoutput = runMonteCarloFluorescence(FMCinput);

%% Post-processing


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

% thetaBeam and phiBeam define the direction of beam center axis, only used if beamtypeflag ~= 1:
% Given in terms of the spherical coordinates theta and phi measured in radians, using the ISO
% convention illustrated at https://en.wikipedia.org/wiki/Spherical_coordinate_system
% Keep in mind that the z-axis in the volumetric plots is shown pointing down, so you want to satisfy 0<=theta<pi/2.
% Examples: theta = 0, phi = 0 means a beam going straight down (positive z direction)
%           theta = pi/4, phi = 0 means a beam going halfway between the positive x and positive z directions.
%           theta = pi/4, phi = -pi/2 means a beam going halfway between the negative y and positive z directions

% A "light collector" can be either an objective lens or a fiber tip.

% theta and phi define the direction that the light collector is facing, defined in the same way as the beam direction using
% ISO spherical coordinates.
% For above the surface (negative z), you'll want to satisfy 0<=theta<=pi/2.
% For below the surface (positive z), such as measuring in transmission, you'll want pi/2<=theta<=pi
% If theta = 0 or pi, phi serves only to rotate the view. If you want the light collector X axis 
% to coincide with the cuboid x axis, use phi = pi/2.

% For time-resolved detection, the results are stored in a number of time
% bins. Specify here the start time, end time and number of bins.
% If you specify nTimeBins = 0, detection is not time-resolved and the light
% collector image output is simply a 2D res by res matrix if the light
% collector is a lens, or a scalar if the light collector is a fiber.
% If you specify nTimeBins > 0, detection is time-resolved and the light
% collector "image" output is a 3D res by res by (nTimeBins+2) matrix if
% the light collector is a lens and a 1D (nTimeBins+2) array if the light
% collector is a fiber. The first time bin stores all photons arriving
% before tStart, the last time bin stores all photons arriving after tEnd,
% and the 2:end-1 time bins equidistantly store the photons arriving
% between tStart and tEnd.

function M = GeometryDefinition_StandardTissue(X,Y,Z,parameters)
tissuedepth = parameters{1};
M = ones(size(X)); % Air
M(Z > tissuedepth) = 3; % "Standard" tissue
end

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

function M = GeometryDefinition_FluorescingCylinder(X,Y,Z,parameters)
cylinderradius  = 0.0100;
M = 17*ones(size(X)); % fill background with fluorescence absorber
M(Y.^2 + (Z - 3*cylinderradius).^2 < cylinderradius^2) = 16; % fluorescer
end

function M = GeometryDefinition_HairExample(X,Y,Z,parameters)
zsurf = 0.02;  % position of gel/skin surface[cm]
epd_thick = 0.01; % thickness of the epidermis [cm]
hair_radius = 0.0075/2; % diameter varies from 17 - 180 micrometers, should increase with colouring and age
hair_bulb_semiminor = 1.7*hair_radius; % [cm]
hair_bulb_semimajor = sqrt(2)*hair_bulb_semiminor;
hair_depth = 0.1; % varies from 0.06-0.3cm
papilla_semiminor = hair_bulb_semiminor*5/12;
papilla_semimajor = sqrt(2)*papilla_semiminor;

M = 2*ones(size(X)); % water (gel)
M(Z > zsurf) = 4; % epidermis
M(Z > zsurf+epd_thick) = 5; % dermis
M(X.^2 + Y.^2 < hair_radius^2 & Z < zsurf+hair_depth) = 10; % hair
M((X/hair_bulb_semiminor).^2 + (Y/hair_bulb_semiminor).^2 + ((Z-(zsurf+hair_depth))/hair_bulb_semimajor).^2 < 1) = 10; % hair
M((X/papilla_semiminor).^2 + (Y/papilla_semiminor).^2 + ((Z-(zsurf+hair_depth+hair_bulb_semimajor-papilla_semimajor))/papilla_semimajor).^2 < 1) = 5; % dermis (papilla)
end

function M = GeometryDefinition_SolderPatchExample(X,Y,Z,parameters)
patch_radius        = 0.218;   	% [cm], cylinder radius
patch_zi_start      = 1;
patch_zi_end        = 5;
vessel_radius       = 0.19;   	% [cm], cylinder radius
water_radius        = 0.15;   	% [cm], cylinder radius
fibre_radius        = 0.04;   	% [cm], cylinder radius

M = ones(size(X)); % fill background with air
M(X.^2 + Y.^2 < patch_radius^2 & Z >= patch_zi_start & Z <= patch_zi_end) = 12; % patch
M(X.^2 + Y.^2 < vessel_radius^2) = 7; % vessel
M(X.^2 + Y.^2 < water_radius^2) = 2; % water
M(X.^2 + Y.^2 < fibre_radius^2) = 11; % fibre
end

function M = GeometryDefinition_ImagingExample(X,Y,Z,parameters)
[nx,ny,~] = size(X);
M = ones(size(X)); % Air background
M(1:(nx*(ny+1)+1):end) = 18; % Set xyz diagonal positions to testscatterer
M(1:(nx*(ny+1)):end) = 18; % Set yz diagonal positions to testscatterer
end

function M = GeometryDefinition_RefractionReflectionExample(X,Y,Z,parameters)
M = ones(size(X)); % Air background
M(Z>0.03) = 2; % Water
M(Z>0.09) = 20; % Reflector
end

function M = GeometryDefinition_stlExample(X,Y,Z,parameters)
xvec = squeeze(X(:,1,1));
yvec = squeeze(Y(1,:,1));
zvec = squeeze(Z(1,1,:));
STLvoxels = VOXELISE(xvec,yvec,zvec,'Dragon_Head.stl');

M = ones(size(X));
M(STLvoxels) = 3;
end
