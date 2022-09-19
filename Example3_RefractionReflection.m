%% Description
% This example illustrates simulation of refraction and reflection
% according to Fresnel's equations. This requires the matchedInterfaces flag
% to be set to false. The geometry consists of three layers: At the bottom
% is a reflector such as metal (n = infinity), at the top is air (n = 1)
% and in between is water (n = 1.3). The light source is an isotropically
% emitting (equally emitting in all directions) source located inside the
% water layer. The rays can be seen to be reflected from the bottom
% interface and also to an extent from the top interface, although some
% light is also refracted out of the water.
%
% This example also shown the optional "useAllCPUs" flag that can be set on
% the MC simulations to slightly increase speed. Default is false, which
% means the solver leaves one CPU unused, making it easier to perform other
% work on the PC while simulations are running. Note that multithreading is
% anyway only supported on Windows. On Mac, useAllCPUs is ignored.
%
% Additionally, the optional "nExamplePaths" parameter is demonstrated, a
% value specifying the number of photons whose paths should be stored and
% shown as lines in a 3D plot after completion.

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

model.G.nx                = 100; % Number of bins in the x direction
model.G.ny                = 100; % Number of bins in the y direction
model.G.nz                = 100; % Number of bins in the z direction
model.G.Lx                = .1; % [cm] x size of simulation cuboid
model.G.Ly                = .1; % [cm] y size of simulation cuboid
model.G.Lz                = .1; % [cm] z size of simulation cuboid

model.G.mediaPropertiesFunc = @mediaPropertiesFunc; % Media properties defined as a function at the end of this file
model.G.geomFunc          = @geometryDefinition; % The distribution of media in the cuboid, also defined as a function at the end of this file.

model = plot(model,'G');

%% Monte Carlo simulation
model.MC.useAllCPUs               = true; % If false, MCmatlab will leave one processor unused. Useful for doing other work on the PC while simulations are running.
model.MC.simulationTimeRequested  = .1; % [min] Time duration of the simulation
model.MC.nExamplePaths            = 100; % (Default: 0) This number of photons will have their paths stored and shown after completion, for illustrative purposes

model.MC.matchedInterfaces        = false; % If false, uses the refractive indices as defined in mediaPropertiesFunc at the end of this file
model.MC.boundaryType             = 1; % 0: No escaping boundaries, 1: All cuboid boundaries are escaping, 2: Top cuboid boundary only is escaping, 3: Top and bottom boundaries are escaping, while the side boundaries are cyclic
model.MC.wavelength               = 532; % [nm] Excitation wavelength, used for determination of optical properties for excitation light

model.MC.lightSource.sourceType   = 1; % 0: Pencil beam, 1: Isotropically emitting line or point source, 2: Infinite plane wave, 3: Laguerre-Gaussian LG01 beam, 4: Radial-factorizable beam (e.g., a Gaussian beam), 5: X/Y factorizable beam (e.g., a rectangular LED emitter)
model.MC.lightSource.xFocus       = 0; % [cm] x position of focus
model.MC.lightSource.yFocus       = 0; % [cm] y position of focus
model.MC.lightSource.zFocus       = 0.04; % [cm] z position of focus

% Execution, do not modify the next line:
model = runMonteCarlo(model);

model = plot(model,'MC');

%% Geometry function(s)
% A geometry function takes as input X,Y,Z matrices as returned by the
% "ndgrid" MATLAB function as well as any parameters the user may have
% provided in the definition of Ginput. It returns the media matrix M,
% containing numerical values indicating the media type (as defined in
% mediaPropertiesFunc) at each voxel location.
function M = geometryDefinition(X,Y,Z,parameters)
  M = ones(size(X)); % Air background
  M(Z>0.03) = 2; % Water
  M(Z>0.09) = 3; % Reflector
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

  j=3;
  mediaProperties(j).name  = 'reflector';
  mediaProperties(j).mua   = 1; % [cm^-1]
  mediaProperties(j).mus   = 1; % [cm^-1]
  mediaProperties(j).g     = 0;
  mediaProperties(j).n     = Inf;
end