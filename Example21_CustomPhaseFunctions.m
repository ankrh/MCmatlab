%% Description
% In this example, we showcase the use of custom phase functions (i.e.,
% other phase functions than Henyey-Greenstein).
%
% We will use the farField output array to see the angular distribution of
% scattered light, so we set model.MC.farFieldRes = 21. Our light source
% will be a pencil beam, going straight down the center column of voxels.
%
% To show the scattering distribution, we use a special geometry in which
% we first fill the cuboid with air, then fill the center column with a
% medium with custom scattering phase function. We choose small values for
% Lx and Ly to ensure that the scattering column is very narrow, so as to
% reduce the likelihood of multiple scattering. To ensure that
% non-scattered photons do not get registered in the farField array we set
% depositionCriteria.minScatterings = 1. See example 19 for how deposition
% criteria work.
%
% For the custom phase function for the scatterer, instead of setting the
% "g" field in the media properties struct we set instead the
% customPhaseFunc field, which has be be a char array that writes out a
% MATLAB-evaluateable expression that is a function of theta. The
% expression does not need to be normalized. MCmatlab will evaluate the
% expression and use an efficient lookup table algorithm to find the
% scattering angle at every scattering event.
%
% In this example, we set the custom phase function to that of Rayleigh
% scattering of unpolarized light. This distribution has equal scattering
% per solid angle in the forward and backward directions, and half of that
% in the sideways directions. This is also what we see in the far field
% plot (figure 9).

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

model.G.nx                = 101; % Number of bins in the x direction
model.G.ny                = 101; % Number of bins in the y direction
model.G.nz                = 150; % Number of bins in the z direction
model.G.Lx                = .02; % [cm] x size of simulation cuboid
model.G.Ly                = .02; % [cm] y size of simulation cuboid
model.G.Lz                = .15; % [cm] z size of simulation cuboid

model.G.mediaPropertiesFunc = @mediaPropertiesFunc; % Media properties defined as a function at the end of this file
model.G.geomFunc          = @geometryDefinition; % Function to use for defining the distribution of media in the cuboid. Defined at the end of this m file.

model = plot(model,'G');

%% Monte Carlo simulation
model.MC.simulationTimeRequested  = .1; % [min] Time duration of the simulation
model.MC.matchedInterfaces        = true; % Assumes all refractive indices are the same
model.MC.boundaryType             = 1; % 0: No escaping boundaries, 1: All cuboid boundaries are escaping, 2: Top cuboid boundary only is escaping, 3: Top and bottom boundaries are escaping, while the side boundaries are cyclic
model.MC.wavelength               = 532; % [nm] Excitation wavelength, used for determination of optical properties for excitation light

model.MC.farFieldRes = 21;

model.MC.lightSource.sourceType   = 0; % 0: Pencil beam, 1: Isotropically emitting line or point source, 2: Infinite plane wave, 3: Laguerre-Gaussian LG01 beam, 4: Radial-factorizable beam (e.g., a Gaussian beam), 5: X/Y factorizable beam (e.g., a rectangular LED emitter)

model.MC.lightSource.xFocus       = 0; % [cm] x position of focus
model.MC.lightSource.yFocus       = 0; % [cm] y position of focus
model.MC.lightSource.zFocus       = model.G.Lz/2; % [cm] z position of focus

model.MC.lightSource.theta        = 0; % [rad] Polar angle of beam center axis
model.MC.lightSource.phi          = 0; % [rad] Azimuthal angle of beam center axis

model.MC.depositionCriteria.minScatterings = 1;


model = runMonteCarlo(model);
model = plot(model,'MC');

%% Geometry function(s) (see readme for details)
function M = geometryDefinition(X,Y,Z,parameters)
  M = ones(size(X)); % Air
  M(51,51,:) = 2; % Scatterer
end

%% Media Properties function (see readme for details)
function mediaProperties = mediaPropertiesFunc(parameters)
  mediaProperties = MCmatlab.mediumProperties;

  j=1;
  mediaProperties(j).name  = 'Air';
  mediaProperties(j).mua   = 1e-8; % [cm^-1]
  mediaProperties(j).mus   = 1e-8; % [cm^-1]
  mediaProperties(j).g     = 1;

  j=2;
  mediaProperties(j).name  = 'Rayleigh scatterer';
  mediaProperties(j).mua   = 1e-8; % [cm^-1]
  mediaProperties(j).mus   = 10; % [cm^-1]
  mediaProperties(j).customPhaseFunc = @func1;
  function phasefunc = func1(lambda,theta)
    phasefunc = 1 + cos(theta)^2;
  end
end