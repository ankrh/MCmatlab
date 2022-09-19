%% Description
% Here we demonstrate how to define beams using sourceType = 1, 4 or 5. The
% geometry is only air, which has negligible scattering and absorption. The
% farFieldRes input will be further illustrated in example 7.
%
% This example will execute five different simulations. Press any key with
% the command window in focus to start the next simulation.
%
% First we demonstrate an isotropically emitting line source with sourceType
% = 1 with length specified using lightSource.emitterLength and focus in the
% middle of the cuboid. If emitterLength is set to 0, the source would be a
% point. The theta and phi parameters determine what direction the line is
% pointing.
%
% Next, we set the focus at x=y=z=0 and theta = 0 so that the beam goes
% straight down. Then we define various simple and complicated beams using
% sourceType 4 and 5. See other examples for more beams with focus placed
% inside the cuboid and for tilted input beams.
%
% In principle, Gaussian beams can be simulated using either sourceType 4 or
% 5 (because exp(-R^2) = exp(-X^2)*exp(-Y^2)).

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
model.G.nz                = 100; % Number of bins in the z direction
model.G.Lx                = .2; % [cm] x size of simulation cuboid
model.G.Ly                = .2; % [cm] y size of simulation cuboid
model.G.Lz                = .1; % [cm] z size of simulation cuboid

model.G.mediaPropertiesFunc = @mediaPropertiesFunc; % Media properties defined as a function at the end of this file
model.G.geomFunc          = @geometryDefinition; % Function to use for defining the distribution of media in the cuboid. Defined at the end of this m file.

%% Monte Carlo simulations
%% Isotropic line emitter
model.MC.simulationTimeRequested  = .1; % [min] Time duration of the simulation
model.MC.farFieldRes              = 200; % (Default: 0) If nonzero, photons that "escape" will have their energies tracked in a 2D angle distribution (theta,phi) array with theta and phi resolutions equal to this number. An "escaping" photon is one that hits the top cuboid boundary (if boundaryType == 2) or any cuboid boundary (if boundaryType == 1) where the medium has refractive index 1.

model.MC.matchedInterfaces        = true; % Assumes all refractive indices are the same
model.MC.boundaryType             = 1; % 0: No escaping boundaries, 1: All cuboid boundaries are escaping, 2: Top cuboid boundary only is escaping, 3: Top and bottom boundaries are escaping, while the side boundaries are cyclic
model.MC.wavelength               = 532; % [nm] Excitation wavelength, used for determination of optical properties for excitation light

model.MC.lightSource.sourceType   = 1; % 0: Pencil beam, 1: Isotropically emitting line or point source, 2: Infinite plane wave, 3: Laguerre-Gaussian LG01 beam, 4: Radial-factorizable beam (e.g., a Gaussian beam), 5: X/Y factorizable beam (e.g., a rectangular LED emitter)
model.MC.lightSource.emitterLength = 0.05; % [cm] (Optional) Length of the isotropically emitting line source (if 0, the source is a point source)

% For an emitter inside the cuboid (point or line) the "focus" position
% just refers to the center of the emitter:
model.MC.lightSource.xFocus       = 0; % [cm] x position of focus
model.MC.lightSource.yFocus       = 0; % [cm] y position of focus
model.MC.lightSource.zFocus       = 0.05; % [cm] z position of focus

model.MC.lightSource.theta        = pi/4; % [rad] Polar angle of light source center axis
model.MC.lightSource.phi          = pi/2; % [rad] Azimuthal angle of light source center axis

% Execution, do not modify the next line:
model = runMonteCarlo(model);
model = plot(model,'MC');
fprintf('Press enter to continue...\n');pause;

%% Top-hat focal plane, Gaussian angular intensity distribution
model.MC.lightSource.xFocus       = 0; % [cm] x position of focus
model.MC.lightSource.yFocus       = 0; % [cm] y position of focus
model.MC.lightSource.zFocus       = 0; % [cm] z position of focus

model.MC.lightSource.sourceType   = 4; % 0: Pencil beam, 1: Isotropically emitting line or point source, 2: Infinite plane wave, 3: Laguerre-Gaussian LG01 beam, 4: Radial-factorizable beam (e.g., a Gaussian beam), 5: X/Y factorizable beam (e.g., a rectangular LED emitter)

model.MC.lightSource.focalPlaneIntensityDistribution.radialDistr = 0; % Radial focal plane intensity distribution - 0: Top-hat, 1: Gaussian, Array: Custom. Doesn't need to be normalized.
model.MC.lightSource.focalPlaneIntensityDistribution.radialWidth = .01; % [cm] Radial focal plane 1/e^2 radius if top-hat or Gaussian or half-width of the full distribution if custom

model.MC.lightSource.angularIntensityDistribution.radialDistr = 1; % Radial angular intensity distribution - 0: Top-hat, 1: Gaussian, 2: Cosine (Lambertian), Array: Custom. Doesn't need to be normalized.
model.MC.lightSource.angularIntensityDistribution.radialWidth = pi/6; % [rad] Radial angular 1/e^2 half-angle if top-hat or Gaussian or half-angle of the full distribution if custom. For a diffraction limited Gaussian beam, this should be set to model.MC.wavelength*1e-9/(pi*model.MC.lightSource.focalPlaneIntensityDistribution.radialWidth*1e-2))

model.MC.lightSource.theta        = 0; % [rad] Polar angle of beam center axis

% Execution, do not modify the next line:
model = runMonteCarlo(model);
model = plot(model,'MC');
fprintf('Press enter to continue...\n');pause;

%% LED-like emitter: Rectangular focal plane, Lambertian angular intensity distribution
model.MC.lightSource.sourceType   = 5; % 0: Pencil beam, 1: Isotropically emitting line or point source, 2: Infinite plane wave, 3: Laguerre-Gaussian LG01 beam, 4: Radial-factorizable beam (e.g., a Gaussian beam), 5: X/Y factorizable beam (e.g., a rectangular LED emitter)

model.MC.lightSource.focalPlaneIntensityDistribution.XDistr = 0; % X focal plane intensity distribution - 0: Top-hat, 1: Gaussian, Array: Custom. Doesn't need to be normalized.
model.MC.lightSource.focalPlaneIntensityDistribution.XWidth = .02; % [cm] X focal plane 1/e^2 radius if top-hat or Gaussian or half-width of the full distribution if custom

model.MC.lightSource.focalPlaneIntensityDistribution.YDistr = 0; % Y focal plane intensity distribution - 0: Top-hat, 1: Gaussian, Array: Custom. Doesn't need to be normalized.
model.MC.lightSource.focalPlaneIntensityDistribution.YWidth = .01; % [cm] Y focal plane 1/e^2 radius if top-hat or Gaussian or half-width of the full distribution if custom

model.MC.lightSource.angularIntensityDistribution.XDistr = 2; % X angular intensity distribution - 0: Top-hat, 1: Gaussian, 2: Cosine (Lambertian), Array: Custom. Doesn't need to be normalized.

model.MC.lightSource.angularIntensityDistribution.YDistr = 2; % Y angular intensity distribution - 0: Top-hat, 1: Gaussian, 2: Cosine (Lambertian), Array: Custom. Doesn't need to be normalized.

model.MC.lightSource.psi          = pi/4; % [rad] (Default: 0) Axial rotation angle of beam, relevant only for XY distributed beams

% Execution, do not modify the next line:
model = runMonteCarlo(model);
model = plot(model,'MC');
fprintf('Press enter to continue...\n');pause;

%% Custom radial focal plane and angular intensity distribution
model.MC.lightSource.sourceType    = 4; % 0: Pencil beam, 1: Isotropically emitting line or point source, 2: Infinite plane wave, 3: Laguerre-Gaussian LG01 beam, 4: Radial-factorizable beam (e.g., a Gaussian beam), 5: X/Y factorizable beam (e.g., a rectangular LED emitter)

% To demonstrate the custom beam definitions, we will construct some
% artificial intensity distributions with some arbitrary analytical
% expressions. In this case the distribution in the focal plane forms a
% ring, while the angular intensity distribution forms four lobes.
customFocalPlaneDistribution = exp(-(linspace(0,20,1000)-4).^2);
model.MC.lightSource.focalPlaneIntensityDistribution.radialDistr = customFocalPlaneDistribution; % Radial focal plane intensity distribution - 0: Top-hat, 1: Gaussian, Array: Custom. Doesn't need to be normalized.
model.MC.lightSource.focalPlaneIntensityDistribution.radialWidth = .025; % [cm] Radial focal plane 1/e^2 radius if top-hat or Gaussian or half-width of the full distribution if custom

customAngularDistribution = 1+cos(linspace(0,7*pi,1000));
model.MC.lightSource.angularIntensityDistribution.radialDistr = customAngularDistribution; % Radial angular intensity distribution - 0: Top-hat, 1: Gaussian, 2: Cosine (Lambertian), Array: Custom. Doesn't need to be normalized.
model.MC.lightSource.angularIntensityDistribution.radialWidth = pi/4; % [rad] Radial angular 1/e^2 half-angle if top-hat or Gaussian or half-angle of the full distribution if custom. For a diffraction limited Gaussian beam, this should be set to model.MC.wavelength*1e-9/(pi*model.MC.lightSource.focalPlaneIntensityDistribution.radialWidth*1e-2))

% Execution, do not modify the next line:
model = runMonteCarlo(model);
model = plot(model,'MC');
fprintf('Press enter to continue...\n');pause;

%% Custom XY focal plane, top-hat X angular intensity distribution and Gaussian Y angular intensity distribution
model.MC.lightSource.sourceType      = 5; % 0: Pencil beam, 1: Isotropically emitting line or point source, 2: Infinite plane wave, 3: Laguerre-Gaussian LG01 beam, 4: Radial-factorizable beam (e.g., a Gaussian beam), 5: X/Y factorizable beam (e.g., a rectangular LED emitter)

% Here, the arbitrarily chosen distribution for the focal plane in the X
% direction consists of two sine wave lobes while the distribution in Y
% consists of three sine wave lobes. For the angular intensity
% distribution, we choose a top-hat along X and a Gaussian along Y.
customXdistribution = sin(linspace(0,2*pi,1000)).^2;
model.MC.lightSource.focalPlaneIntensityDistribution.XDistr = customXdistribution; % X focal plane intensity distribution - 0: Top-hat, 1: Gaussian, Array: Custom. Doesn't need to be normalized.
model.MC.lightSource.focalPlaneIntensityDistribution.XWidth = .02; % [cm] X focal plane 1/e^2 radius if top-hat or Gaussian or half-width of the full distribution if custom

customYdistribution = sin(linspace(0,3*pi,1000)).^2;
model.MC.lightSource.focalPlaneIntensityDistribution.YDistr = customYdistribution; % Y focal plane intensity distribution - 0: Top-hat, 1: Gaussian, Array: Custom. Doesn't need to be normalized.
model.MC.lightSource.focalPlaneIntensityDistribution.YWidth = .01; % [cm] Y focal plane 1/e^2 radius if top-hat or Gaussian or half-width of the full distribution if custom

model.MC.lightSource.angularIntensityDistribution.XDistr = 0; % X angular intensity distribution - 0: Top-hat, 1: Gaussian, 2: Cosine (Lambertian), Array: Custom. Doesn't need to be normalized.
model.MC.lightSource.angularIntensityDistribution.XWidth = pi/8; % [rad] X angular 1/e^2 half-angle if top-hat or Gaussian or half-angle of the full distribution if custom

model.MC.lightSource.angularIntensityDistribution.YDistr = 1; % Y angular intensity distribution - 0: Top-hat, 1: Gaussian, 2: Cosine (Lambertian), Array: Custom. Doesn't need to be normalized.
model.MC.lightSource.angularIntensityDistribution.YWidth = pi/8; % [rad] Y angular 1/e^2 half-angle if top-hat or Gaussian or half-angle of the full distribution if custom

model.MC.lightSource.psi             = -pi/4; % [rad] (Default: 0) Axial rotation angle of beam, relevant only for XY distributed beams

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
  M = ones(size(X)); % Air
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
end