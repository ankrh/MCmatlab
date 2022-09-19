%% Description
% This example illustrates that deposition criteria can be specified to determine when
% photon packets should contribute to output data, for example one could
% specify the number of scattering events must be at least 1 and at most 3.
%
% The list of criteria are:
% Minimum and maximum # of scattering events
% Minimum and maximum # of refraction events (for non-matched interfaces)
% Minimum and maximum # of reflection events (for non-matched interfaces)
% Minimum and maximum # of interface transition events
% 
% There is an implied AND between all these critera.
%
% When a photon packet fulfills the criteria, its "weight loss"
% (absorption) will be counted in the normalized fluence rate 3D array, in
% the image 2D array, in the 2D boundary irradiance arrays and in the
% angular-resolved far field data. Also it will be registered as a
% potential example photon path.
%
% The default criteria are zeros for the minima and infinities for the
% maxima. That means that by default, all photon packets deposit weight
% into the output arrays regardless of the number of scatterings,
% refractions etc.
%
% Changing deposition criteria does not change the underlying propagation
% and absorption of photons packets in the simulation. It only means that
% data will only be registered (deposited) into the output arrays when the
% photon packets fulfill the criteria.
%
% In this example, we use the same geometry and settings as example 1,
% except that we add some photon paths and we choose only to plot data for
% photons that have been scattered exactly once by setting
% model.MC.depositionCriteria.minScatterings = 1 and
% model.MC.depositionCriteria.maxScatterings = 1. Notice how you can see the linear
% paths of photons scattered away from the center column.
%
% See the readme.md file for all the properties of the depositionCriteria
% object.

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
model.MC.wavelength               = 532; % [nm] Excitation wavelength, used for determination of optical properties for excitation light

model.MC.lightSource.sourceType   = 0; % 0: Pencil beam, 1: Isotropically emitting line or point source, 2: Infinite plane wave, 3: Laguerre-Gaussian LG01 beam, 4: Radial-factorizable beam (e.g., a Gaussian beam), 5: X/Y factorizable beam (e.g., a rectangular LED emitter)

model.MC.lightSource.xFocus       = 0; % [cm] x position of focus
model.MC.lightSource.yFocus       = 0; % [cm] y position of focus
model.MC.lightSource.zFocus       = model.G.Lz/2; % [cm] z position of focus

model.MC.lightSource.theta        = 0; % [rad] Polar angle of beam center axis
model.MC.lightSource.phi          = 0; % [rad] Azimuthal angle of beam center axis

model.MC.nExamplePaths = 25;
model.MC.depositionCriteria.minScatterings = 1; % Don't collect data for non-scattered photons
model.MC.depositionCriteria.maxScatterings = 1; % Don't collect data for photons with more than 1 scattering

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
  zSurface = 0.03;
  M = ones(size(X)); % Air
  M(Z > zSurface) = 2; % "Standard" tissue
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

  j=2;
  mediaProperties(j).name  = 'standard tissue';
  mediaProperties(j).mua   = 1; % [cm^-1]
  mediaProperties(j).mus   = 10; % [cm^-1]
  mediaProperties(j).g     = 0;
end