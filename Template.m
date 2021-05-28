%% Description
%
%

%% Geometry definition
model = MCmatlab.model;

model.G.silentMode        = false; % (Default: false) Disables command window text and progress indication

model.G.nx                = 100; % Number of bins in the x direction
model.G.ny                = 100; % Number of bins in the y direction
model.G.nz                = 100; % Number of bins in the z direction
model.G.Lx                = .1; % [cm] x size of simulation cuboid
model.G.Ly                = .1; % [cm] y size of simulation cuboid
model.G.Lz                = .1; % [cm] z size of simulation cuboid

% model.G.mediaPropParams   = {0.6}; % Cell array containing any additional parameters to be passed to the getMediaProperties function

model.G.mediaPropertiesFunc = @mediaPropertiesFunc; % Media properties defined as a function at the end of this file
model.G.geomFunc          = @geometryDefinition_StandardTissue; % Function to use for defining the distribution of media in the cuboid. Defined at the end of this m file.
% model.G.geomFuncParams    = {0.03}; % Cell array containing any additional parameters to pass into the geometry function, such as media depths, inhomogeneity positions, radii etc.

plot(model,'G');

%% Monte Carlo simulation
% model = reset(model,'MC'); % Only necessary if you want to run this section repeatedly, re-using previous G data
% 
% model.MC.useGPU = false; % (Default: false) Use CUDA acceleration for NVIDIA GPUs
% 
% model.MC.simulationTimeRequested  = .1; % [min] Time duration of the simulation
% model.MC.nPhotonsRequested        = 1e5; % # of photons to launch
% 
% % model.MC.silentMode               = false; % (Default: false) Disables command window text and progress indication
% % model.MC.useAllCPUs               = true; % (Default: false) If false, MCmatlab will leave one processor unused. Useful for doing other work on the PC while simulations are running.
% % model.MC.calcNFR                  = true; % (Default: true) If true, the 3D fluence rate output array NFR will be calculated. Set to false if you have a light collector and you're only interested in the image output.
% % model.MC.calcNFRdet               = false; % (Default: false) If true, the 3D fluence rate output array NFRdet will be calculated. Only photons that end up on the light collector are counted in NFRdet.
% % model.MC.nExamplePaths            = 100; % (Default: 0) This number of photons will have their paths stored and shown after completion, for illustrative purposes
% % model.MC.farFieldRes              = 50; % (Default: 0) If nonzero, photons that "escape" will have their energies tracked in a 2D angle distribution (theta,phi) array with theta and phi resolutions equal to this number. An "escaping" photon is one that hits the top cuboid boundary (if boundaryType == 2) or any cuboid boundary (if boundaryType == 1) where the medium has refractive index 1.
% 
% model.MC.matchedInterfaces        = true; % (Default: true) If false, uses the refractive indices as defined in mediaPropertiesFunc at the end of this file
% model.MC.smoothingLengthScale     = model.G.Lx*2; % [cm] The characteristic length scale to smoothe the map of surface normal vectors over, for use when simulating refraction and reflection angles
% model.MC.boundaryType             = 1; % 0: No escaping boundaries, 1: All cuboid boundaries are escaping, 2: Top cuboid boundary only is escaping
% model.MC.wavelength               = 450; % [nm] Excitation wavelength, used for determination of optical properties for excitation light
% 
% model.MC.beam.beamType            = 5; % 0: Pencil beam, 1: Isotropically emitting line or point source, 2: Infinite plane wave, 3: Laguerre-Gaussian LG01 beam, 4: Radial-factorizable beam (e.g., a Gaussian beam), 5: X/Y factorizable beam (e.g., a rectangular LED emitter)
% model.MC.beam.emitterLength       = 0.05; % [cm] (Optional) Length of the isotropically emitting line source (if 0, the source is a point source)
% model.MC.beam.NF.radialDistr      = exp(-(linspace(0,20,1000)-4).^2); % Radial near field distribution - 0: Top-hat, 1: Gaussian, Array: Custom. Doesn't need to be normalized.
% model.MC.beam.NF.radialWidth      = .025; % [cm] Radial near field 1/e^2 radius if top-hat or Gaussian or half-width of the full distribution if custom
% model.MC.beam.NF.XDistr           = sin(linspace(0,2*pi,1000)).^2; % X near field distribution - 0: Top-hat, 1: Gaussian, Array: Custom. Doesn't need to be normalized.
% model.MC.beam.NF.XWidth           = .02; % [cm] X near field 1/e^2 radius if top-hat or Gaussian or half-width of the full distribution if custom
% model.MC.beam.NF.YDistr           = sin(linspace(0,3*pi,1000)).^2; % Y near field distribution - 0: Top-hat, 1: Gaussian, Array: Custom. Doesn't need to be normalized.
% model.MC.beam.NF.YWidth           = .01; % [cm] Y near field 1/e^2 radius if top-hat or Gaussian or half-width of the full distribution if custom
% model.MC.beam.FF.radialDistr      = 1+cos(linspace(0,7*pi,1000)); % Radial far field distribution - 0: Top-hat, 1: Gaussian, 2: Cosine (Lambertian), Array: Custom. Doesn't need to be normalized.
% model.MC.beam.FF.radialWidth      = pi/4; % [rad] Radial far field 1/e^2 half-angle if top-hat or Gaussian or half-angle of the full distribution if custom. For a diffraction limited Gaussian beam, this should be set to model.MC.wavelength*1e-9/(pi*model.MC.beam.NF.radialWidth*1e-2))
% model.MC.beam.FF.XDistr           = 0; % X far field distribution - 0: Top-hat, 1: Gaussian, 2: Cosine (Lambertian), Array: Custom. Doesn't need to be normalized.
% model.MC.beam.FF.XWidth           = pi/8; % [rad] X far field 1/e^2 half-angle if top-hat or Gaussian or half-angle of the full distribution if custom
% model.MC.beam.FF.YDistr           = 1; % Y far field distribution - 0: Top-hat, 1: Gaussian, 2: Cosine (Lambertian), Array: Custom. Doesn't need to be normalized.
% model.MC.beam.FF.YWidth           = pi/8; % [rad] Y far field 1/e^2 half-angle if top-hat or Gaussian or half-angle of the full distribution if custom
%
% model.MC.beam.xFocus              = 0; % [cm] x position of focus
% model.MC.beam.yFocus              = 0; % [cm] y position of focus
% model.MC.beam.zFocus              = 0; % [cm] z position of focus
%
% model.MC.beam.psi                 = -pi/4; % [rad] (Default: 0) Axial rotation angle of beam, relevant only for XY distributed beams
%
% model.MC.beam.theta               = 0; % [rad] Polar angle of beam center axis
% model.MC.beam.phi                 = 0; % [rad] Azimuthal angle of beam center axis
% 
% % model.MC.P                        = 4; % [W] Incident pulse peak power (in case of infinite plane waves, only the power incident upon the cuboid's top surface)
% % model.MC.FRinitial = zeros(model.G.nx,model.G.ny,model.G.nz); % [W/cm^2] Initial guess for the intensity distribution, to be used for fluence rate dependent simulations
% % model.MC.FRdepIterations = 20;
% 
% % model.MC.useLightCollector      = true;
% % model.MC.LC.x         = 0; % [cm] x position of either the center of the objective lens focal plane or the fiber tip
% % model.MC.LC.y         = 0; % [cm] y position
% % model.MC.LC.z         = 0.03; % [cm] z position
% % 
% % model.MC.LC.theta     = 0; % [rad] Polar angle of direction the light collector is facing
% % model.MC.LC.phi       = pi/2; % [rad] Azimuthal angle of direction the light collector is facing
% % 
% % model.MC.LC.f         = .2; % [cm] Focal length of the objective lens (if light collector is a fiber, set this to Inf).
% % model.MC.LC.diam      = .1; % [cm] Diameter of the light collector aperture. For an ideal thin lens, this is 2*f*tan(asin(NA)).
% % model.MC.LC.fieldSize = .1; % [cm] Field Size of the imaging system (diameter of area in object plane that gets imaged). Only used for finite f.
% % model.MC.LC.NA        = 0.22; % [-] Fiber NA. Only used for infinite f.
% % 
% % model.MC.LC.res       = 50; % X and Y resolution of light collector in pixels, only used for finite f
% % 
% % % model.MC.LC.tStart    = -1e-13; % [s] Start of the detection time-of-flight interval
% % % model.MC.LC.tEnd      = 5e-12; % [s] End of the detection time-of-flight interval
% % % model.MC.LC.nTimeBins = 30; % (Default: 0) Number of bins between tStart and tEnd. If zero, the measurement is not time-resolved.
% 
% % Execution, do not modify the next line:
% model = runMonteCarlo(model);
% 
% plot(model,'MC');

%% Fluorescence Monte Carlo
% model = reset(model,'FMC'); % Only necessary if you want to run this section repeatedly, re-using previous G and MC data
% 
% model.FMC.useGPU = false; % (Default: false) Use CUDA acceleration for NVIDIA GPUs
% 
% model.FMC.simulationTimeRequested  = .1; % [min] Time duration of the simulation
% model.FMC.nPhotonsRequested        = 1e5; % # of photons to launch
% 
% % model.FMC.silentMode             = false; % (Default: false) Disables command window text and progress indication
% % model.FMC.useAllCPUs             = true; % (Default: false) If false, MCmatlab will leave one processor unused. Useful for doing other work on the PC while simulations are running.
% % model.FMC.calcNFR                = true; % (Default: true) If true, the 3D fluence rate output array NFR will be calculated. Set to false if you have a light collector and you're only interested in the image output.
% % model.FMC.calcNFRdet             = false; % (Default: false) If true, the 3D fluence rate output array NFRdet will be calculated. Only photons that end up on the light collector are counted in NFRdet.
% % model.FMC.nExamplePaths          = 100; % (Default: 0) This number of photons will have their paths stored and shown after completion, for illustrative purposes
% % model.FMC.farFieldRes            = 50; % (Default: 0) If nonzero, photons that "escape" will have their energies tracked in a 2D angle distribution (theta,phi) array with theta and phi resolutions equal to this number. An "escaping" photon is one that hits the top cuboid boundary (if boundaryType == 2) or any cuboid boundary (if boundaryType == 1) where the medium has refractive index 1.
% 
% model.FMC.matchedInterfaces        = true; % (Default: true) If false, uses the refractive indices as defined in mediaPropertiesFunc at the end of this file
% model.FMC.smoothingLengthScale     = model.G.Lx*2; % [cm] The characteristic length scale to smoothe the map of surface normal vectors over, for use when simulating refraction and reflection angles
% model.FMC.boundaryType             = 1; % 0: No escaping boundaries, 1: All cuboid boundaries are escaping, 2: Top cuboid boundary only is escaping
% model.FMC.wavelength               = 450; % [nm] Excitation wavelength, used for determination of optical properties for excitation light
% 
% % model.FMC.useLightCollector      = true;
% % model.FMC.LC.x         = 0; % [cm] x position of either the center of the objective lens focal plane or the fiber tip
% % model.FMC.LC.y         = 0; % [cm] y position
% % model.FMC.LC.z         = 0.03; % [cm] z position
% % 
% % model.FMC.LC.theta     = 0; % [rad] Polar angle of direction the light collector is facing
% % model.FMC.LC.phi       = pi/2; % [rad] Azimuthal angle of direction the light collector is facing
% % 
% % model.FMC.LC.f         = .2; % [cm] Focal length of the objective lens (if light collector is a fiber, set this to Inf).
% % model.FMC.LC.diam      = .1; % [cm] Diameter of the light collector aperture. For an ideal thin lens, this is 2*f*tan(asin(NA)).
% % model.FMC.LC.fieldSize = .1; % [cm] Field Size of the imaging system (diameter of area in object plane that gets imaged). Only used for finite f.
% % model.FMC.LC.NA        = 0.22; % [-] Fiber NA. Only used for infinite f.
% % 
% % model.FMC.LC.res       = 50; % X and Y resolution of light collector in pixels, only used for finite f
% 
% % Execution, do not modify the next line:
% model = runMonteCarlo(model,'fluorescence');
% 
% plot(model,'FMC');

%% Heat simulation
% model = reset(model,'HS'); % Only necessary if you want to run this section repeatedly, re-using previous G, MC and/or FMC data
% 
% model.HS.useGPU = false; % (Default: false) Use CUDA acceleration for NVIDIA GPUs
% % model.HS.silentMode          = false; % (Default: false) Disables command window text and progress indication
% % model.HS.useAllCPUs          = true; % (Default: false) If false, MCmatlab will leave one processor unused. Useful for doing other work on the PC while simulations are running.
% % model.HS.makeMovie           = true; % (Default: false) Requires silentMode = false.
% % model.HS.largeTimeSteps      = true; % (Default: false) If true, calculations will be faster, but some voxel temperatures may be slightly less precise. Test for yourself whether this precision is acceptable for your application.
% % model.HS.deferMovieWrite = false;
% 
% model.HS.heatBoundaryType    = 0; % 0: Insulating boundaries, 1: Constant-temperature boundaries (heat-sinked)
% model.HS.durationOn          = 0.001; % [s] Pulse on-duration
% model.HS.durationOff         = 0.004; % [s] Pulse off-duration
% model.HS.durationEnd         = 0.02; % [s] Non-illuminated relaxation time to add to the end of the simulation to let temperature diffuse after the pulse train
% model.HS.Tinitial            = 37; % [deg C] Initial temperature, can be scalar or 3D array
% 
% model.HS.nPulses             = 5; % Number of consecutive pulses, each with an illumination phase and a diffusion phase. If simulating only illumination or only diffusion, use n_pulses = 1.
% 
% model.HS.plotTempLimits      = [37 100]; % [deg C] Expected range of temperatures, used only for setting the color scale in the plot
% model.HS.nUpdates            = 100; % Number of times data is extracted for plots during each pulse. A minimum of 1 update is performed in each phase (2 for each pulse consisting of an illumination phase and a diffusion phase)
% % model.HS.mediaPropRecalcPeriod = 5; % Every N updates, the media properties will be recalculated (including, if needed, re-running MC and FMC steps)
% 
% model.HS.slicePositions      = [.5 0.6 1]; % (Default: [0.5 1 1]) Relative slice positions [x y z] for the 3D plots on a scale from 0 to 1
% model.HS.tempSensorPositions = [0 0 0.038
%                               0 0 0.04
%                               0 0 0.042
%                               0 0 0.044]; % (Default: []) Each row is a temperature sensor's absolute [x y z] coordinates. Leave the matrix empty ([]) to disable temperature sensors.
% 
% % Execution, do not modify the next line:
% model = simulateHeatDistribution(model);
% 
% plot(model,'HS');

%% Post-processing

%% Geometry function(s)
% A geometry function takes as input X,Y,Z matrices as returned by the
% "ndgrid" MATLAB function as well as any parameters the user may have
% provided in the definition of model.G. It returns the media matrix M,
% containing numerical values indicating the media type (as defined in
% mediaPropertiesFunc) at each voxel location.
function M = geometryDefinition_StandardTissue(X,Y,Z,parameters)
tissuedepth = 0.03;
M = ones(size(X)); % Air
M(Z > tissuedepth) = 3; % "Standard" tissue
end

function M = geometryDefinition_BloodVessel(X,Y,Z,parameters)
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

function M = geometryDefinition_FluorescingCylinder(X,Y,Z,parameters)
cylinderradius  = 0.0100;
M = 17*ones(size(X)); % fill background with fluorescence absorber
M(Y.^2 + (Z - 3*cylinderradius).^2 < cylinderradius^2) = 16; % fluorescer
end

function M = geometryDefinition_HairExample(X,Y,Z,parameters)
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

function M = geometryDefinition_SolderPatchExample(X,Y,Z,parameters)
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

function M = geometryDefinition_TimeTaggingExample(X,Y,Z,parameters)
[nx,ny,~] = size(X);
M = ones(size(X)); % Air background
M(1:(nx*(ny+1)+1):end) = 18; % Set xyz diagonal positions to testscatterer
M(1:(nx*(ny+1)):end) = 18; % Set yz diagonal positions to testscatterer
end

function M = geometryDefinition_RefractionReflectionExample(X,Y,Z,parameters)
M = ones(size(X)); % Air background
M(Z>0.03) = 2; % Water
M(Z>0.09) = 20; % Reflector
end

%% Media Properties function
% The media properties function defines all the optical and thermal
% properties of the media involved by constructing and returning a
% "mediaProperties" struct with various fields. As its input, the function
% takes the wavelength as well as any other parameters you might specify
% above in the model file, for example parameters that you might loop over
% in a for loop. Dependence on excitation fluence rate FR, temperature T or
% fractional heat damage FD can be specified as in examples 12-15.
function mediaProperties = mediaPropertiesFunc(wavelength,parameters)
j=1;
mediaProperties(j).name  = 'air';
mediaProperties(j).mua   = 1e-8;
mediaProperties(j).mus   = 1e-8;
mediaProperties(j).g     = 1;
mediaProperties(j).n     = 1;
mediaProperties(j).VHC   = 1.2e-3;
mediaProperties(j).TC    = 0; % Real value is 2.6e-4, but we set it to zero to neglect the heat transport to air

j=2;
mediaProperties(j).name  = 'standard tissue';
mediaProperties(j).mua   = 1;
mediaProperties(j).mus   = 100;
mediaProperties(j).g     = 0.9;
mediaProperties(j).n     = 1.3;
mediaProperties(j).VHC   = 3391*1.109e-3;
mediaProperties(j).TC    = 0.37e-2;
end