%% Description
% In this example, we show the use of multiple light collectors.

%% MCmatlab abbreviations
% G: Geometry, MC: Monte Carlo, FMC: Fluorescence Monte Carlo, HS: Heat
% simulation, M: Media array, FR: Fluence rate, FD: Fractional damage.
%
% There are also some optional abbreviations you can use when referencing
% object/variable names: LS = lightSource, Dets = detectors, FPID =
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
model.MC.simulationTimeRequested  = .01; % [min] Time duration of the simulation
model.MC.nExamplePaths            = 100; % (Default: 0) This number of photons will have their paths stored and shown after completion, for illustrative purposes

model.MC.matchedInterfaces        = true; % Assumes all refractive indices are the same
model.MC.boundaryType             = 1; % 0: No escaping boundaries, 1: All cuboid boundaries are escaping, 2: Top cuboid boundary only is escaping, 3: Top and bottom boundaries are escaping, while the side boundaries are cyclic
model.MC.wavelength               = 450; % [nm] Excitation wavelength, used for determination of optical properties for excitation light

model.MC.lightSource.sourceType   = 2; % 0: Pencil beam, 1: Isotropically emitting line or point source, 2: Infinite plane wave, 3: Laguerre-Gaussian LG01 beam, 4: Radial-factorizable beam (e.g., a Gaussian beam), 5: X/Y factorizable beam (e.g., a rectangular LED emitter)
model.MC.lightSource.theta        = 0; % [rad] Polar angle of beam center axis
model.MC.lightSource.phi          = 0; % [rad] Azimuthal angle of beam center axis

model.MC.detectors(1).x      = 0; % [cm] x position of either the center of the objective lens focal plane or the fiber tip
model.MC.detectors(1).y      = 0; % [cm] y position
model.MC.detectors(1).z      = -0.03; % [cm] z position
model.MC.detectors(1).theta  = 0;
model.MC.detectors(1).phi    = 0;
model.MC.detectors(1).psi    = 0;
model.MC.detectors(1).Xsize  = .03; % [cm] Size of the detectors in the X direction.
model.MC.detectors(1).Ysize  = .06; % [cm] Size of the detectors in the X direction.
model.MC.detectors(1).shape  = 'Ellipse'; % Can be either ellipse or rectangle.
model.MC.detectors(1).res    = 3;

model.MC.detectors(2).x      = 0.02; % [cm] x position of either the center of the objective lens focal plane or the fiber tip
model.MC.detectors(2).y      = 0; % [cm] y position
model.MC.detectors(2).z      = -0.03; % [cm] z position
model.MC.detectors(2).theta  = 0;
model.MC.detectors(2).phi    = 0;
model.MC.detectors(2).psi    = 0;
model.MC.detectors(2).Xsize  = .06; % [cm] Size of the detectors in the X direction.
model.MC.detectors(2).Ysize  = .02; % [cm] Size of the detectors in the X direction.
model.MC.detectors(2).shape  = 'Rectangle'; % Can be either ellipse or rectangle.

model = runMonteCarlo(model);
model = plot(model,'MC');

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
