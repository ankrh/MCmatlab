%% Description
% This example shows the use of the helper functions MCmatlab.plotAzFz()
% and MCmatlab.plotRadialReflectanceTransmittance(). They generate plots
% similar to those made by MCML, available at omlc.org. Use the log scale
% button to see the details in the wings of the distributions.
% 
% The radially resolved reflectance and transmittance plots only make sense
% if your problem is radially symmetric around x = y = 0. Furthermore, when
% launching simulations for these plots, it's advised to use an odd number
% for nx and ny, to ensure a data point at a radial distance of 0.

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
model.G.Lx                = 5; % [cm] x size of simulation cuboid
model.G.Ly                = 5; % [cm] y size of simulation cuboid
model.G.Lz                = 1; % [cm] z size of simulation cuboid

model.G.mediaPropertiesFunc = @mediaPropertiesFunc; % Media properties defined as a function at the end of this file
model.G.geomFunc          = @geometryDefinition; % Function to use for defining the distribution of media in the cuboid. Defined at the end of this m file.

model = plot(model,'G');

%% Monte Carlo simulation
model.MC.simulationTimeRequested  = .1; % [min] Time duration of the simulation
model.MC.matchedInterfaces        = true; % Assumes all refractive indices are the same
model.MC.boundaryType             = 1; % 0: No escaping boundaries, 1: All cuboid boundaries are escaping, 2: Top cuboid boundary only is escaping, 3: Top and bottom boundaries are escaping, while the side boundaries are cyclic
model.MC.wavelength               = 532; % [nm] Excitation wavelength, used for determination of optical properties for excitation light

model.MC.lightSource.sourceType   = 0; % 0: Pencil beam, 1: Isotropically emitting line or point source, 2: Infinite plane wave, 3: Laguerre-Gaussian LG01 beam, 4: Radial-factorizable beam (e.g., a Gaussian beam), 5: X/Y factorizable beam (e.g., a rectangular LED emitter)

model.MC.lightSource.xFocus       = 0; % [cm] x position of focus
model.MC.lightSource.yFocus       = 0; % [cm] y position of focus
model.MC.lightSource.zFocus       = model.G.Lz/2; % [cm] z position of focus
model.MC.lightSource.theta        = 0; % [rad] Polar angle of beam center axis
model.MC.lightSource.phi          = 0; % [rad] Azimuthal angle of beam center axis

model = runMonteCarlo(model);
model = plot(model,'MC');

% Make MCML-like Az, Fz, reflectance and transmission plots
MCmatlab.plotAzFz(model);
MCmatlab.plotRadialReflectanceTransmittance(model);

%% Geometry function(s) (see readme for details)
function M = geometryDefinition(X,Y,Z,parameters)
  M = ones(size(X)); % Tissue 1
  M(Z > 0.2) = 2; % Tissue 2
end

%% Media Properties function (see readme for details)
function mediaProperties = mediaPropertiesFunc(parameters)
  mediaProperties = MCmatlab.mediumProperties;

  j=1;
  mediaProperties(j).name  = 'Tissue 1';
  mediaProperties(j).mua   = 0.4; % Absorption coefficient [cm^-1]
  mediaProperties(j).mus   = 100; % Scattering coefficient [cm^-1]
  mediaProperties(j).g     = 0.9; % Henyey-Greenstein scattering anisotropy

  j=2;
  mediaProperties(j).name  = 'Tissue 2';
  mediaProperties(j).mua   = 0.2; % Absorption coefficient [cm^-1]
  mediaProperties(j).mus   = 100; % Scattering coefficient [cm^-1]
  mediaProperties(j).g     = 0.9; % Henyey-Greenstein scattering anisotropy
end