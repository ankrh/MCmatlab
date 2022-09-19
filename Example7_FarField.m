%% Description
% Here, we demonstrate calculation and plotting of the far field angular
% distribution of the "escaped" photons, both for excitation light and for
% fluorescence light. It is enabled by specifying the optional
% "farFieldRes" parameter, which serves also to specify the resolution you
% want of the far field distribution. The geometry is almost the same
% fluorescing cylinder as example 5, but now illuminated by a tilted
% Gaussian beam and now with only one type of fluorescer.
%
% In the far field of the excitation light, you can see that they primarily
% escape in downward-pointing directions (in transmission), while the far
% field distribution of fluorescence light indicates that fluorescence
% mostly escapes on the ends of the cylinder, with a small amount of light
% coming out of the cylinder in the minus z direction (upwards).
%
% Note that MCmatlab distinguishes between "escaped" photons and "killed"
% photons. An "escaping" photon is one that hits the top cuboid boundary
% (if boundaryType == 2) or any cuboid boundary (if boundaryType == 1)
% where the medium has refractive index 1. A "killed" photon is one that
% strays too far from the main cuboid (5 times further than the cuboid
% dimensions).
%
% If boundaryType == 1, there are no "killed" photons since no photons can
% travel outside the cuboid, and the fraction of light absorbed in the
% cuboid plus the fraction of light escaping equals 1.

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
model.G.nz                = 101; % Number of bins in the z direction
model.G.Lx                = .1; % [cm] x size of simulation cuboid
model.G.Ly                = .1; % [cm] y size of simulation cuboid
model.G.Lz                = .1; % [cm] z size of simulation cuboid

model.G.mediaPropertiesFunc = @mediaPropertiesFunc; % Media properties defined as a function at the end of this file
model.G.geomFunc          = @geometryDefinition; % Function to use for defining the distribution of media in the cuboid. Defined at the end of this m file.

model = plot(model,'G');

%% Monte Carlo simulation
model.MC.simulationTimeRequested  = .1; % [min] Time duration of the simulation
model.MC.farFieldRes              = 50; % (Default: 0) If nonzero, photons that "escape" will have their energies tracked in a 2D angle distribution (theta,phi) array with theta and phi resolutions equal to this number. An "escaping" photon is one that hits the top cuboid boundary (if boundaryType == 2) or any cuboid boundary (if boundaryType == 1) where the medium has refractive index 1.

model.MC.matchedInterfaces        = true; % Assumes all refractive indices are the same
model.MC.boundaryType             = 1; % 0: No escaping boundaries, 1: All cuboid boundaries are escaping, 2: Top cuboid boundary only is escaping, 3: Top and bottom boundaries are escaping, while the side boundaries are cyclic
model.MC.wavelength               = 450; % [nm] Excitation wavelength, used for determination of optical properties for excitation light

model.MC.lightSource.sourceType   = 4; % 0: Pencil beam, 1: Isotropically emitting line or point source, 2: Infinite plane wave, 3: Laguerre-Gaussian LG01 beam, 4: Radial-factorizable beam (e.g., a Gaussian beam), 5: X/Y factorizable beam (e.g., a rectangular LED emitter)
model.MC.lightSource.focalPlaneIntensityDistribution.radialDistr = 1; % Radial focal plane intensity distribution - 0: Top-hat, 1: Gaussian, Array: Custom. Doesn't need to be normalized.
model.MC.lightSource.focalPlaneIntensityDistribution.radialWidth = .015; % [cm] Radial focal plane 1/e^2 radius if top-hat or Gaussian or half-width of the full distribution if custom
model.MC.lightSource.angularIntensityDistribution.radialDistr = 1; % Radial angular intensity distribution - 0: Top-hat, 1: Gaussian, 2: Cosine (Lambertian), Array: Custom. Doesn't need to be normalized.
model.MC.lightSource.angularIntensityDistribution.radialWidth = 15/180*pi; % [rad] Radial angular 1/e^2 half-angle if top-hat or Gaussian or half-angle of the full distribution if custom. For a diffraction limited Gaussian beam, this should be set to model.MC.wavelength*1e-9/(pi*model.MC.lightSource.focalPlaneIntensityDistribution.radialWidth*1e-2))
model.MC.lightSource.xFocus       = 0; % [cm] x position of focus
model.MC.lightSource.yFocus       = 0; % [cm] y position of focus
model.MC.lightSource.zFocus       = 0.03; % [cm] z position of focus
model.MC.lightSource.theta        = pi/6; % [rad] Polar angle of beam center axis
model.MC.lightSource.phi          = -pi/2; % [rad] Azimuthal angle of beam center axis

% Execution, do not modify the next line:
model = runMonteCarlo(model);

model = plot(model,'MC');

%% Fluorescence Monte Carlo
model.FMC.simulationTimeRequested  = 0.1; % [min] Time duration of the simulation
model.FMC.farFieldRes              = 50; % (Default: 0) If nonzero, photons that "escape" will have their energies tracked in a 2D angle distribution (theta,phi) array with theta and phi resolutions equal to this number. An "escaping" photon is one that hits the top cuboid boundary (if boundaryType == 2) or any cuboid boundary (if boundaryType == 1) where the medium has refractive index 1.

model.FMC.matchedInterfaces        = true; % Assumes all refractive indices are the same
model.FMC.boundaryType             = 1; % 0: No escaping boundaries, 1: All cuboid boundaries are escaping, 2: Top cuboid boundary only is escaping, 3: Top and bottom boundaries are escaping, while the side boundaries are cyclic
model.FMC.wavelength               = 550; % [nm] Fluorescence wavelength, used for determination of optical properties for fluorescence light

% Execution, do not modify the next line:
model = runMonteCarlo(model,'fluorescence');

model = plot(model,'FMC');

%% Geometry function(s)
% A geometry function takes as input X,Y,Z matrices as returned by the
% "ndgrid" MATLAB function as well as any parameters the user may have
% provided in the definition of model.G. It returns the media matrix M,
% containing numerical values indicating the media type (as defined in
% mediaPropertiesFunc) at each voxel location.
function M = geometryDefinition(X,Y,Z,parameters)
  cylinderradius  = 0.0100;
  M = 1*ones(size(X)); % fill background with fluorescence absorber
  M(Y.^2 + (Z - 3*cylinderradius).^2 < cylinderradius^2) = 2; % fluorescer
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
  mediaProperties(j).name  = 'fluorescence absorber';
  if(wavelength<500)
    mediaProperties(j).mua = 1; % [cm^-1]
    mediaProperties(j).mus = 100; % [cm^-1]
    mediaProperties(j).g   = 0.9;
  else
    mediaProperties(j).mua = 100; % [cm^-1]
    mediaProperties(j).mus = 100; % [cm^-1]
    mediaProperties(j).g   = 0.9;
  end

  j=2;
  mediaProperties(j).name  = 'fluorescer';
  if(wavelength<500)
    mediaProperties(j).mua = 100; % [cm^-1]
    mediaProperties(j).mus = 100; % [cm^-1]
    mediaProperties(j).g   = 0.9;
  else
    mediaProperties(j).mua = 1; % [cm^-1]
    mediaProperties(j).mus = 100; % [cm^-1]
    mediaProperties(j).g   = 0.9;
  end

  % Only one of PY and QY may be defined:
  mediaProperties(j).PY   = 0.5; % Fluorescence power yield (ratio of power emitted to power absorbed)
  % mediaProperties(j).QY   = 0.6; % Fluorescence quantum yield (ratio of photons emitted to photons absorbed)
end
