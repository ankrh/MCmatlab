%% Description
% This example illustrates refraction and reflection using the smoothed and
% interpolated Sobel filter method published by A. Tran and S. Jacques in
% J. of Biomedical Optics, 25(2), 025001 (2020)
% https://doi.org/10.1117/1.JBO.25.2.025001).
%
% This example reproduces figure 3 of that paper. A hemispherical droplet
% of water is surrounded by air and an infinite plane wave is incident on
% it. The two media, of course, have different refractive indices and the
% normal vectors of the droplet surface are calculated using a Sobel filter
% and smoothed at a length scale of model.MC.smoothingLengthScale.
%
% The focusing effect of the droplet is clearly visible when looking at the
% photon paths. The probability of reflection is greater near the edges due
% to the grazing angles of incidence. You will see some voxelization
% artifacts in the light distribution near the edges of the droplet, which
% are unavoidable with this method but will be less significant if the
% resolution is increased.

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
model = MCmatlab.model;

model.G.nx                = 201; % Number of bins in the x direction
model.G.ny                = 201; % Number of bins in the y direction
model.G.nz                = 101; % Number of bins in the z direction
model.G.Lx                = .1; % [cm] x size of simulation cuboid
model.G.Ly                = .1; % [cm] y size of simulation cuboid
model.G.Lz                = .05; % [cm] z size of simulation cuboid

model.G.mediaPropertiesFunc = @mediaPropertiesFunc; % Media properties defined as a function at the end of this file
model.G.geomFunc          = @geometryDefinition; % The distribution of media in the cuboid, also defined as a function at the end of this file.

plot(model,'G');

%% Monte Carlo simulation
model.MC.useAllCPUs               = true; % If false, MCmatlab will leave one processor unused. Useful for doing other work on the PC while simulations are running.

model.MC.simulationTimeRequested  = .1; % [min] Time duration of the simulation
model.MC.nExamplePaths            = 100; % (Default: 0) This number of photons will have their paths stored and shown after completion, for illustrative purposes

model.MC.matchedInterfaces        = false; % If false, uses the refractive indices as defined in mediaPropertiesFunc at the end of this file
model.MC.smoothingLengthScale     = model.G.Lx*10; % [cm] The characteristic length scale to smoothe the map of surface normal vectors over, for use when simulating refraction and reflection angles

model.MC.boundaryType             = 1; % 0: No escaping boundaries, 1: All cuboid boundaries are escaping, 2: Top cuboid boundary only is escaping, 3: Top and bottom boundaries are escaping, while the side boundaries are cyclic
model.MC.wavelength               = 532; % [nm] Excitation wavelength, used for determination of optical properties for excitation light

model.MC.lightSource.sourceType   = 2; % 0: Pencil beam, 1: Isotropically emitting line or point source, 2: Infinite plane wave, 3: Laguerre-Gaussian LG01 beam, 4: Radial-factorizable beam (e.g., a Gaussian beam), 5: X/Y factorizable beam (e.g., a rectangular LED emitter)

model.MC.lightSource.focalPlaneIntensityDistribution.radialDistr = 1; % Radial focal plane intensity distribution - 0: Top-hat, 1: Gaussian, Array: Custom. Doesn't need to be normalized.
model.MC.lightSource.focalPlaneIntensityDistribution.radialWidth = .01; % [cm] Radial focal plane 1/e^2 radius if top-hat or Gaussian or half-width of the full distribution if custom
model.MC.lightSource.angularIntensityDistribution.radialDistr = 1; % Radial angular intensity distribution - 0: Top-hat, 1: Gaussian, 2: Cosine (Lambertian), Array: Custom. Doesn't need to be normalized.
model.MC.lightSource.angularIntensityDistribution.radialWidth = 0; % [rad] Radial angular 1/e^2 half-angle if top-hat or Gaussian or half-angle of the full distribution if custom. For a diffraction limited Gaussian beam, this should be set to model.MC.wavelength*1e-9/(pi*model.MC.lightSource.focalPlaneIntensityDistribution.radialWidth*1e-2))
model.MC.lightSource.xFocus       = 0; % [cm] x position of focus
model.MC.lightSource.yFocus       = 0; % [cm] y position of focus
model.MC.lightSource.zFocus       = model.G.Lz; % [cm] z position of focus
model.MC.lightSource.theta        = 0; % [rad] Polar angle of beam center axis
model.MC.lightSource.phi          = 0; % [rad] Azimuthal angle of beam center axis

% Execution, do not modify the next line:
model = runMonteCarlo(model);

plot(model,'MC');

%% Geometry function(s)
% A geometry function takes as input X,Y,Z matrices as returned by the
% "ndgrid" MATLAB function as well as any parameters the user may have
% provided in the definition of Ginput. It returns the media matrix M,
% containing numerical values indicating the media type (as defined in
% mediaPropertiesFunc) at each voxel location.
function M = geometryDefinition(X,Y,Z,parameters)
  M = ones(size(X)); % Air background
  M(X.^2 + Y.^2 + (Z - Z(end)).^2 < 0.04^2) = 2; % Water
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
  mediaProperties(j).mua   = 1e-8; % [cm^-1]
  mediaProperties(j).mus   = 1e-8; % [cm^-1]
  mediaProperties(j).g     = 1;
  mediaProperties(j).n     = 1;

  j=2;
  mediaProperties(j).name  = 'water';
  mediaProperties(j).mua   = 0.00036; % [cm^-1]
  mediaProperties(j).mus   = 10; % [cm^-1]
  mediaProperties(j).g     = 1.0;
  mediaProperties(j).n     = 1.3;
end