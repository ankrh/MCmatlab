%% Description
% In this example, we show that MCmatlab can launch excitation photons from
% a user-specified spatial and spectral distribution of sources, similar to
% how fluorescence photons are launched. The distribution may be wavelength
% dependent. It is stored in model.MC.sourceDistribution and must be 3- or
% 4-dimensional (nx, ny, nz) or (nx, ny, nz, nLambda).

% The distribution of emitters gets plotted in figure 7.

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

model.G.nx                = 101; % Number of bins in the x direction
model.G.ny                = 101; % Number of bins in the y direction
model.G.nz                = 150; % Number of bins in the z direction
model.G.Lx                = .1; % [cm] x size of simulation cuboid
model.G.Ly                = .1; % [cm] y size of simulation cuboid
model.G.Lz                = .15; % [cm] z size of simulation cuboid

model.G.mediaPropertiesFunc = @mediaPropertiesFunc; % Media properties defined as a function at the end of this file
model.G.geomFunc          = @geometryDefinition; % Function to use for defining the distribution of media in the cuboid. Defined at the end of this m file.

model = plot(model,'G');

%% Monte Carlo simulation
model.MC.simulationTimeRequested  = .1; % [min] Time duration of the simulation
model.MC.matchedInterfaces        = true; % Assumes all refractive indices are the same
model.MC.boundaryType             = 1; % 0: No escaping boundaries, 1: All cuboid boundaries are escaping, 2: Top cuboid boundary only is escaping, 3: Top and bottom boundaries are escaping, while the side boundaries are cyclic
model.MC.wavelength               = linspace(400,600,41); % [nm] Excitation wavelength, used for determination of optical properties for excitation light

% For our spatial-spectral 4D array for our source distribution, we
% arbitrarily choose a sphere buried in the tissue, with a Gaussian
% spectral distribution. The distribution will automatically be normalized
% once runMonteCarlo() is called:
[X,Y,Z,lambda] = ndgrid(model.G.x, model.G.y, model.G.z, model.MC.wavelength);
model.MC.sourceDistribution = exp(-(lambda-500).^2/50^2).*(X.^2+Y.^2+(Z-0.1).^2 < 0.01^2);


model = runMonteCarlo(model);
model = plot(model,'MC');

%% Geometry function(s) (see readme for details)
function M = geometryDefinition(X,Y,Z,parameters)
  zSurface = 0.03;
  M = ones(size(X)); % Air
  M(Z > zSurface) = 2; % "Standard" tissue
end

%% Media Properties function (see readme for details)
function mediaProperties = mediaPropertiesFunc(parameters)
  mediaProperties = MCmatlab.mediumProperties;

  j=1;
  mediaProperties(j).name  = 'air';
  mediaProperties(j).mua   = 1e-8; % Absorption coefficient [cm^-1]
  mediaProperties(j).mus   = 1e-8; % Scattering coefficient [cm^-1]
  mediaProperties(j).g     = 1; % Henyey-Greenstein scattering anisotropy

  j=2;
  mediaProperties(j).name  = 'standard tissue';
  mediaProperties(j).mua   = 1; % Absorption coefficient [cm^-1]
  mediaProperties(j).mus   = 100; % Scattering coefficient [cm^-1]
  mediaProperties(j).g     = 0.9; % Henyey-Greenstein scattering anisotropy
end