%% Description
% In this example, we show the use of multiple detectors. Detectors with
% res > 1 are spatially resolved and plots will be generated in figures
% 101, 102, 103 etc. with the irradiance distributions.

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

model.MC.simulationTimeRequested  = .1; % [min] Time duration of the simulation

model.MC.matchedInterfaces        = true; % Assumes all refractive indices are the same
model.MC.boundaryType             = 1; % 0: No escaping boundaries, 1: All cuboid boundaries are escaping, 2: Top cuboid boundary only is escaping, 3: Top and bottom boundaries are escaping, while the side boundaries are cyclic
model.MC.wavelength               = 450; % [nm] Excitation wavelength, used for determination of optical properties for excitation light

model.MC.lightSource.sourceType   = 2; % 0: Pencil beam, 1: Isotropically emitting line or point source, 2: Infinite plane wave, 3: Laguerre-Gaussian LG01 beam, 4: Radial-factorizable beam (e.g., a Gaussian beam), 5: X/Y factorizable beam (e.g., a rectangular LED emitter)
model.MC.lightSource.theta        = 0; % [rad] Polar angle of beam center axis
model.MC.lightSource.phi          = 0; % [rad] Azimuthal angle of beam center axis

j = 1;
model.MC.detectors(j).x      = -0.051; % [cm] x position of the center of the detector
model.MC.detectors(j).y      = 0.01; % [cm] y position
model.MC.detectors(j).z      = 0.05; % [cm] z position
model.MC.detectors(j).theta  = pi/2;
model.MC.detectors(j).phi    = 0;
model.MC.detectors(j).psi    = pi/2;
model.MC.detectors(j).Xsize  = .15; % [cm] Size of the detectors in the X direction.
model.MC.detectors(j).Ysize  = .10; % [cm] Size of the detectors in the X direction.
model.MC.detectors(j).shape  = 'Ellipse'; % Can be either circle, ellipse, square or rectangle. Circle is an alias for ellipse and square is an alias for rectangle.
model.MC.detectors(j).res    = 30;

j = 2;
model.MC.detectors(j).x      = 0.05; % [cm] x position of the center of the detector
model.MC.detectors(j).y      = 0.05; % [cm] y position
model.MC.detectors(j).z      = -0.005; % [cm] z position
model.MC.detectors(j).theta  = 0;
model.MC.detectors(j).phi    = 0;
model.MC.detectors(j).psi    = pi/2;
model.MC.detectors(j).Xsize  = .06; % [cm] Size of the detectors in the X direction.
model.MC.detectors(j).Ysize  = .02; % [cm] Size of the detectors in the X direction.
model.MC.detectors(j).shape  = 'Rectangle'; % Can be either circle, ellipse, square or rectangle. Circle is an alias for ellipse and square is an alias for rectangle.
model.MC.detectors(j).res    = 30;

j = 3;
model.MC.detectors(j).x      = 0; % [cm] x position of the center of the detector
model.MC.detectors(j).y      = 0; % [cm] y position
model.MC.detectors(j).z      = 0.11; % [cm] z position
model.MC.detectors(j).theta  = pi;
model.MC.detectors(j).phi    = 0;
model.MC.detectors(j).psi    = 0;
model.MC.detectors(j).Xsize  = .15; % [cm] Size of the detectors in the X direction.
model.MC.detectors(j).Ysize  = .15; % [cm] Size of the detectors in the X direction.
model.MC.detectors(j).shape  = 'Rectangle'; % Can be either circle, ellipse, square or rectangle. Circle is an alias for ellipse and square is an alias for rectangle.
model.MC.detectors(j).res    = 1;

model = runMonteCarlo(model);

model = plot(model,'MC');

%% Geometry function(s) (see readme for details)
function M = geometryDefinition(X,Y,Z,parameters)
  M = ones(size(X));
end

%% Media Properties function (see readme for details)
function mediaProperties = mediaPropertiesFunc(parameters)
  mediaProperties = MCmatlab.mediumProperties;

  j=1;
  mediaProperties(j).name  = 'standard tissue';
  mediaProperties(j).mua   = 1; % Absorption coefficient [cm^-1]
  mediaProperties(j).mus   = 100; % Scattering coefficient [cm^-1]
  mediaProperties(j).g     = 0.9; % Henyey-Greenstein scattering anisotropy
end
