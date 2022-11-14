%% Description
% In this example, simulation of fluorescence (luminescence) is shown. The
% test geometry is a fluorescing cylinder in which excitation light is
% predominantly absorbed embedded in a block of medium in which
% fluorescence light is predominantly absorbed. The geometry is illuminated
% with an infinite plane wave.
%
% This example also shows detection of the light exiting the cuboid,
% separately for excitation light and for fluorescence light. Although most
% of the fluorescence light is absorbed in the medium surrounding the
% cylinder, some of it escapes to the detector, showing a slightly blurred
% image of the cylinder.
%
% To use a light collector, the cuboid boundary type towards the detector
% has to be set to "escaping". Additionally, the voxels touching that
% boundary must have a refractive index of 1.
%
% The "nExamplePaths" parameter (see example 3) is used for both the
% excitation and fluorescence simulations, showing paths of both kinds of
% photons.

%% MCmatlab abbreviations
% G: Geometry, MC: Monte Carlo, FMC: Fluorescence Monte Carlo, HS: Heat
% simulation, M: Media array, FR: Fluence rate, FD: Fractional damage.
%
% There are also some optional abbreviations you can use when referencing
% object/variable names: LS = lightSource, LC = lightCollector, FPID =
% focalPlaneIntensityDistribution, AID = angularIntensityDistribution, NI =
% normalizedIrradiance, NFR = normalizedFluenceRate.
%
% For example, "model.MC.LS.FPID.radialDistr" is the same as
% "model.MC.lightSource.focalPlaneIntensityDistribution.radialDistr"

%% Geometry definition
MCmatlab.closeMCmatlabFigures();
model = MCmatlab.model;

model.G.nx                = 100; % Number of bins in the x direction
model.G.ny                = 100; % Number of bins in the y direction
model.G.nz                = 100; % Number of bins in the z direction
model.G.Lx                = .1; % [cm] x size of simulation cuboid
model.G.Ly                = .1; % [cm] y size of simulation cuboid
model.G.Lz                = .1; % [cm] z size of simulation cuboid

model.G.mediaPropertiesFunc = @mediaPropertiesFunc; % Media properties defined as a function at the end of this file
model.G.geomFunc          = @geometryDefinition; % Function to use for defining the distribution of media in the cuboid. Defined at the end of this m file.

model = plot(model,'G');

%% Monte Carlo simulation
model.MC.useAllCPUs               = true; % If false, MCmatlab will leave one processor unused. Useful for doing other work on the PC while simulations are running.
model.MC.simulationTimeRequested  = .1; % [min] Time duration of the simulation
model.MC.nExamplePaths            = 100; % (Default: 0) This number of photons will have their paths stored and shown after completion, for illustrative purposes

model.MC.matchedInterfaces        = true; % Assumes all refractive indices are the same
model.MC.boundaryType             = 1; % 0: No escaping boundaries, 1: All cuboid boundaries are escaping, 2: Top cuboid boundary only is escaping, 3: Top and bottom boundaries are escaping, while the side boundaries are cyclic
model.MC.wavelength               = 450; % [nm] Excitation wavelength, used for determination of optical properties for excitation light

model.MC.lightSource.sourceType   = 2; % 0: Pencil beam, 1: Isotropically emitting line or point source, 2: Infinite plane wave, 3: Laguerre-Gaussian LG01 beam, 4: Radial-factorizable beam (e.g., a Gaussian beam), 5: X/Y factorizable beam (e.g., a rectangular LED emitter)
model.MC.lightSource.theta        = 0; % [rad] Polar angle of beam center axis
model.MC.lightSource.phi          = 0; % [rad] Azimuthal angle of beam center axis

model.MC.useLightCollector        = true;

model.MC.lightCollector.x         = 0; % [cm] x position of either the center of the objective lens focal plane or the fiber tip
model.MC.lightCollector.y         = 0; % [cm] y position
model.MC.lightCollector.z         = 0.03; % [cm] z position

model.MC.lightCollector.theta     = 0; % [rad] Polar angle of direction the light collector is facing
model.MC.lightCollector.phi       = pi/2; % [rad] Azimuthal angle of direction the light collector is facing

model.MC.lightCollector.f         = .2; % [cm] Focal length of the objective lens (if light collector is a fiber, set this to Inf).
model.MC.lightCollector.diam      = .1; % [cm] Diameter of the light collector aperture. For an ideal thin lens, this is 2*f*tan(asin(NA)).
model.MC.lightCollector.fieldSize = .1; % [cm] Field Size of the imaging system (diameter of area in object plane that gets imaged). Only used for finite f.
model.MC.lightCollector.NA        = 0.22; % [-] Fiber NA. Only used for infinite f.

model.MC.lightCollector.res       = 50; % X and Y resolution of light collector in pixels, only used for finite f


model = runMonteCarlo(model);
model = plot(model,'MC');

%% Fluorescence Monte Carlo
model.FMC.useAllCPUs              = true; % If false, MCmatlab will leave one processor unused. Useful for doing other work on the PC while simulations are running.
model.FMC.simulationTimeRequested = .1; % [min] Time duration of the simulation
model.FMC.nExamplePaths           = 100; % (Default: 0) This number of photons will have their paths stored and shown after completion, for illustrative purposes

model.FMC.matchedInterfaces       = true; % Assumes all refractive indices are the same
model.FMC.boundaryType            = 1; % 0: No escaping boundaries, 1: All cuboid boundaries are escaping, 2: Top cuboid boundary only is escaping, 3: Top and bottom boundaries are escaping, while the side boundaries are cyclic
model.FMC.wavelength              = 900; % [nm] Fluorescence wavelength, used for determination of optical properties for fluorescence light

model.FMC.useLightCollector       = true;

model.FMC.lightCollector.x        = 0; % [cm] x position of either the center of the objective lens focal plane or the fiber tip
model.FMC.lightCollector.y        = 0; % [cm] y position
model.FMC.lightCollector.z        = 0.03; % [cm] z position

model.FMC.lightCollector.theta    = 0; % [rad] Polar angle of direction the light collector is facing
model.FMC.lightCollector.phi      = pi/2; % [rad] Azimuthal angle of direction the light collector is facing

model.FMC.lightCollector.f        = .2; % [cm] Focal length of the objective lens (if light collector is a fiber, set this to Inf).
model.FMC.lightCollector.diam     = .1; % [cm] Diameter of the light collector aperture. For an ideal thin lens, this is 2*f*tan(asin(NA)).
model.FMC.lightCollector.fieldSize = .1; % [cm] Field Size of the imaging system (diameter of area in object plane that gets imaged). Only used for finite f.
model.FMC.lightCollector.NA       = 0.22; % [-] Fiber NA. Only used for infinite f.

model.FMC.lightCollector.res      = 50; % X and Y resolution of light collector in pixels, only used for finite f


model = runMonteCarlo(model,'fluorescence');
model = plot(model,'FMC');

%% Geometry function(s) (see readme for details)
function M = geometryDefinition(X,Y,Z,parameters)
  cylinderradius  = 0.0100;
  M = ones(size(X)); % fill background with fluorescence absorber
  M(Y.^2 + (Z - 3*cylinderradius).^2 < cylinderradius^2) = 2; % fluorescer
end

%% Media Properties function (see readme for details)
function mediaProperties = mediaPropertiesFunc(parameters)
  mediaProperties = MCmatlab.mediumProperties;

  j=1;
  mediaProperties(j).name  = 'fluorescence absorber';
  mediaProperties(j).mua = @func1; % [cm^-1]
  function mua = func1(wavelength)
    if(wavelength<500)
      mua = 1; % [cm^-1]
    else
      mua = 100; % [cm^-1]
    end
  end
  mediaProperties(j).mus = 100; % [cm^-1]
  mediaProperties(j).g   = 0.9;

  j=2;
  mediaProperties(j).name  = 'fluorescer';
  mediaProperties(j).mua = @func2; % [cm^-1]
  function mua = func2(wavelength)
    if(wavelength<500)
      mua = 100; % [cm^-1]
    else
      mua = 1; % [cm^-1]
    end
  end
  mediaProperties(j).mus = 100; % [cm^-1]
  mediaProperties(j).g   = 0.9;

  mediaProperties(j).QY   = 0.4; % Fluorescence quantum yield
end
