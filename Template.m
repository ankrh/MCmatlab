addpath([fileparts(matlab.desktop.editor.getActiveFilename) '/helperfuncs']); % The helperfuncs folder is added to the path for the duration of this MATLAB session

%% Geometry definition
clear Ginput
Ginput.silentMode        = false; % (Default: false) Disables command window text and progress indication
Ginput.matchedInterfaces = true; % (Default: true) If true, assumes all refractive indices are 1. If false, uses the refractive indices defined in getMediaProperties
Ginput.boundaryType      = 1; % 0: No escaping boundaries, 1: All cuboid boundaries are escaping, 2: Top cuboid boundary only is escaping

Ginput.wavelength        = 450; % [nm] Excitation wavelength, used for determination of optical properties for excitation light
% Ginput.wavelength_f      = 550; % [nm] Fluorescence wavelength, used for determination of optical properties for fluorescence light
% Ginput.mediaPropParams   = {0.6}; % Cell array containing any additional parameters to be passed to the getMediaProperties function

Ginput.nx                = 100; % Number of bins in the x direction
Ginput.ny                = 100; % Number of bins in the y direction
Ginput.nz                = 100; % Number of bins in the z direction
Ginput.Lx                = .1; % [cm] x size of simulation cuboid
Ginput.Ly                = .1; % [cm] y size of simulation cuboid
Ginput.Lz                = .1; % [cm] z size of simulation cuboid

Ginput.GeomFunc          = @GeometryDefinition_FluorescingCylinder; % Function to use for defining the distribution of media in the cuboid. Defined at the end of this m file.
% Ginput.GeomFuncParams    = {0.03}; % Cell array containing any additional parameters to pass into the geometry function, such as media depths, inhomogeneity positions, radii etc.

% Execution, do not modify the next two lines:
Goutput = defineGeometry(Ginput);
plotMCmatlabGeom(Goutput);

%% Monte Carlo simulation
% clear MCinput
% MCinput.simulationTime           = .1; % [min] Time duration of the simulation
% MCinput.nPhotons                 = 1e5; % # of photons to launch
% 
% % MCinput.silentMode               = false; % (Default: false) Disables command window text and progress indication
% % MCinput.useAllCPUs               = true; % (Default: false) If false, MCmatlab will leave one processor unused. Useful for doing other work on the PC while simulations are running.
% % MCinput.calcF                    = true; % (Default: true) If true, the 3D fluence rate output matrix F will be calculated. Set to false if you have a light collector and you're only interested in the Image output.
% % MCinput.calcFdet                 = false; % (Default: false) If true, the 3D fluence rate output matrix Fdet will be calculated. Only photons that end up on the light collector are counted in Fdet.
% % MCinput.nExamplePaths            = 100; % (Default: 0) This number of photons will have their paths stored and shown after completion, for illustrative purposes
% % MCinput.farfieldRes              = 50; % (Default: 0) If nonzero, photons that "escape" will have their energies tracked in a 2D angle distribution (theta,phi) array with theta and phi resolutions equal to this number. An "escaping" photon is one that hits the top cuboid boundary (if boundaryType == 2) or any cuboid boundary (if boundaryType == 1) where the medium has refractive index 1.
% 
% MCinput.Beam.beamType            = 2; % 0: Pencil beam, 1: Isotropically emitting point source, 2: Infinite plane wave, 3: Gaussian focus, Gaussian far field beam, 4: Gaussian focus, top-hat far field beam, 5: Top-hat focus, Gaussian far field beam, 6: Top-hat focus, top-hat far field beam, 7: Laguerre-Gaussian LG01 beam
% MCinput.Beam.nearFieldType       = 2; % 0: Gaussian, 1: Circular top-hat, 2: Square top-hat
% MCinput.Beam.farFieldType        = 2; % 0: Gaussian, 1: Circular top-hat, 2: Cosine distribution (Lambertian)
% MCinput.Beam.xFocus              = 0; % [cm] x position of focus
% MCinput.Beam.yFocus              = 0; % [cm] y position of focus
% MCinput.Beam.zFocus              = Ginput.Lz/2; % [cm] z position of focus
% MCinput.Beam.theta               = 0; % [rad] Polar angle of beam center axis
% MCinput.Beam.phi                 = 0; % [rad] Azimuthal angle of beam center axis
% MCinput.Beam.waist               = 0.005; % [cm] Beam waist 1/e^2 radius
% MCinput.Beam.divergence          = 5/180*pi; % [rad] Beam divergence 1/e^2 half-angle of beam (for a diffraction limited Gaussian beam, this is G.wavelength*1e-9/(pi*MCinput.Beam.waist*1e-2))
% 
% % MCinput.LightCollector.x         = 0; % [cm] x position of either the center of the objective lens focal plane or the fiber tip
% % MCinput.LightCollector.y         = 0; % [cm] y position
% % MCinput.LightCollector.z         = 0.03; % [cm] z position
% % 
% % MCinput.LightCollector.theta     = 0; % [rad] Polar angle of direction the light collector is facing
% % MCinput.LightCollector.phi       = pi/2; % [rad] Azimuthal angle of direction the light collector is facing
% % 
% % MCinput.LightCollector.f         = .2; % [cm] Focal length of the objective lens (if light collector is a fiber, set this to Inf).
% % MCinput.LightCollector.diam      = .1; % [cm] Diameter of the light collector aperture. For an ideal thin lens, this is 2*f*tan(asin(NA)).
% % MCinput.LightCollector.FieldSize = .1; % [cm] Field Size of the imaging system (diameter of area in object plane that gets imaged). Only used for finite f.
% % MCinput.LightCollector.NA        = 0.22; % [-] Fiber NA. Only used for infinite f.
% % 
% % MCinput.LightCollector.res       = 50; % X and Y resolution of light collector in pixels, only used for finite f
% % 
% % % MCinput.LightCollector.tStart    = -1e-13; % [s] Start of the detection time interval
% % % MCinput.LightCollector.tEnd      = 5e-12; % [s] End of the detection time interval
% % % MCinput.LightCollector.nTimeBins = 30; % (Default: 0) Number of bins between tStart and tEnd. If zero, the measurement is not time-resolved.
% 
% % Execution, do not modify the next three lines:
% MCinput.G = Goutput;
% MCoutput = runMonteCarlo(MCinput);
% plotMCmatlab(MCinput,MCoutput);

%% Fluorescence Monte Carlo
% clear FMCinput
% FMCinput.simulationTime           = .1; % [min] Time duration of the simulation
% FMCinput.nPhotons                 = 1e5; % # of photons to launch
% 
% % FMCinput.silentMode               = false; % (Default: false) Disables command window text and progress indication
% % FMCinput.useAllCPUs               = true; % (Default: false) If false, MCmatlab will leave one processor unused. Useful for doing other work on the PC while simulations are running.
% % FMCinput.calcF                    = true; % (Default: true) If true, the 3D fluence rate output matrix F will be calculated. Set to false if you have a light collector and you're only interested in the Image output.
% % FMCinput.calcFdet                 = false; % (Default: false) If true, the 3D fluence rate output matrix Fdet will be calculated. Only photons that end up on the light collector are counted in Fdet.
% % FMCinput.nExamplePaths            = 100; % (Default: 0) This number of photons will have their paths stored and shown after completion, for illustrative purposes
% % FMCinput.farfieldRes              = 50; % (Default: 0) If nonzero, photons that "escape" will have their energies tracked in a 2D angle distribution (theta,phi) array with theta and phi resolutions equal to this number. An "escaping" photon is one that hits the top cuboid boundary (if boundaryType == 2) or any cuboid boundary (if boundaryType == 1) where the medium has refractive index 1.
% 
% % FMCinput.LightCollector.x         = 0; % [cm] x position of either the center of the objective lens focal plane or the fiber tip
% % FMCinput.LightCollector.y         = 0; % [cm] y position
% % FMCinput.LightCollector.z         = 0.03; % [cm] z position
% % 
% % FMCinput.LightCollector.theta     = 0; % [rad] Polar angle of direction the light collector is facing
% % FMCinput.LightCollector.phi       = pi/2; % [rad] Azimuthal angle of direction the light collector is facing
% % 
% % FMCinput.LightCollector.f         = .2; % [cm] Focal length of the objective lens (if light collector is a fiber, set this to Inf).
% % FMCinput.LightCollector.diam      = .1; % [cm] Diameter of the light collector aperture. For an ideal thin lens, this is 2*f*tan(asin(NA)).
% % FMCinput.LightCollector.FieldSize = .1; % [cm] Field Size of the imaging system (diameter of area in object plane that gets imaged). Only used for finite f.
% % FMCinput.LightCollector.NA        = 0.22; % [-] Fiber NA. Only used for infinite f.
% % 
% % FMCinput.LightCollector.res       = 50; % X and Y resolution of light collector in pixels, only used for finite f
% 
% % Execution, do not modify the next four lines:
% FMCinput.G = Goutput;
% FMCinput.MCoutput = MCoutput;
% FMCoutput = runMonteCarloFluorescence(FMCinput);
% plotMCmatlabFluorescence(FMCinput,FMCoutput);

%% Heat simulation
% % HSinput.silentMode          = false; % (Default: false) Disables command window text and progress indication
% % HSinput.useAllCPUs          = true; % (Default: false) If false, MCmatlab will leave one processor unused. Useful for doing other work on the PC while simulations are running.
% % HSinput.makeMovie           = true; % (Default: false) Requires silentMode = false.
% % HSinput.largeTimeSteps      = true; % (Default: false) If true, calculations will be faster, but some voxel temperatures may be slightly less precise. Test for yourself whether this precision is acceptable for your application.
% 
% HSinput.heatBoundaryType    = 0; % 0: Insulating boundaries, 1: Constant-temperature boundaries (heat-sinked)
% HSinput.P                   = 4; % [W] Incident pulse peak power (in case of infinite plane waves, only the power incident upon the cuboid's top surface)
% HSinput.durationOn          = 0.001; % [s] Pulse on-duration
% HSinput.durationOff         = 0.004; % [s] Pulse off-duration
% HSinput.durationEnd         = 0.02; % [s] Non-illuminated relaxation time to add to the end of the simulation to let temperature diffuse after the pulse train
% HSinput.initialTemp         = 37; % [deg C] Initial temperature
% 
% HSinput.nPulses             = 5; % Number of consecutive pulses, each with an illumination phase and a diffusion phase. If simulating only illumination or only diffusion, use n_pulses = 1.
% 
% HSinput.plotTempLimits      = [37 100]; % [deg C] Expected range of temperatures, used only for setting the color scale in the plot
% HSinput.nUpdates            = 100; % Number of times data is extracted for plots during each pulse. A minimum of 1 update is performed in each phase (2 for each pulse consisting of an illumination phase and a diffusion phase)
% HSinput.slicePositions      = [.5 0.6 1]; % (Default: [0.5 1 1]) Relative slice positions [x y z] for the 3D plots on a scale from 0 to 1
% HSinput.tempSensorPositions = [0 0 0.038
%                               0 0 0.04
%                               0 0 0.042
%                               0 0 0.044]; % (Default: []) Each row is a temperature sensor's absolute [x y z] coordinates. Leave the matrix empty ([]) to disable temperature sensors.
% 
% % Execution, do not modify the next four lines:
% HSinput.G = Goutput;
% HSinput.MCoutput = MCoutput;
% HSoutput = simulateHeatDistribution(HSinput);
% plotMCmatlabHeat(HSinput,HSoutput);

%% Post-processing

%% Geometry function(s)
% A geometry function takes as input X,Y,Z matrices as returned by the
% "ndgrid" MATLAB function as well as any parameters the user may have
% provided in the definition of Ginput. It returns the media matrix M,
% containing numerical values indicating the media type (as defined in
% getMediaProperties) at each voxel location.
function M = GeometryDefinition_StandardTissue(X,Y,Z,parameters)
tissuedepth = 0.03;
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

function M = GeometryDefinition_TimeTaggingExample(X,Y,Z,parameters)
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

