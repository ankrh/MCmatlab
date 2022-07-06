%% Description
% In this example, we demonstrate that model outputs can be combined with
% each other using the combineModels() function. This can be useful, for
% example, when you want to have multiple input light sources, such as
% beams. It can also be useful for heavy simulations that require a lot of
% statistics, so you might have run the same simulation on multiple
% different PCs and would like to combine the results.

% This example is almost the same as example 1, except that we start by
% naming the model "model1" instead of the usual "model". It has a shifted
% y position of the input pencil beam and a power (model1.MC.P) of 1
% watt.

% Then we copy this model into a new model called "model2", and we give
% model2 a new beam y position and a different power, 0.1 watts.

% We run the Monte Carlo simulation on both of these models individually,
% and then use the combineModels() function to combine the results into one
% model called "combinedModel".

% This combined model object now contains properties where the input objects'...
% - nPhotons and simulationTime have been added
% - examplePaths have been concatenated
% - image, normalizedFluenceRate, farField and all the normalizedIrradiance
%   arrays have been set to weighted means
% If all the light sources (model.MC.LS) and input powers (model.MC.P) are
% identical we assume that the objects are being combined simply to gather
% more statistics, so the weights used for the mean are the numbers of
% photons in each simulation, nPhotons. If they are not identical, we
% assume that we are combining simulations with different light sources, as
% is the case in this example, and the weights used are the powers,
% model.MC.P.

% When plotting the combined model output, you see the two pencil beams in
% the cuboid, with one having five times the power of the other.

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
model1 = MCmatlab.model;

model1.G.nx                = 101; % Number of bins in the x direction
model1.G.ny                = 101; % Number of bins in the y direction
model1.G.nz                = 150; % Number of bins in the z direction
model1.G.Lx                = .1; % [cm] x size of simulation cuboid
model1.G.Ly                = .1; % [cm] y size of simulation cuboid
model1.G.Lz                = .15; % [cm] z size of simulation cuboid

model1.G.mediaPropertiesFunc = @mediaPropertiesFunc; % Media properties defined as a function at the end of this file
model1.G.geomFunc          = @geometryDefinition; % Function to use for defining the distribution of media in the cuboid. Defined at the end of this m file.

plot(model1,'G');

%% Monte Carlo simulation
model1.MC.simulationTimeRequested  = .1; % [min] Time duration of the simulation
model1.MC.matchedInterfaces        = true; % Assumes all refractive indices are the same
model1.MC.boundaryType             = 1; % 0: No escaping boundaries, 1: All cuboid boundaries are escaping, 2: Top cuboid boundary only is escaping, 3: Top and bottom boundaries are escaping, while the side boundaries are cyclic
model1.MC.wavelength               = 532; % [nm] Excitation wavelength, used for determination of optical properties for excitation light

model1.MC.lightSource.sourceType   = 0; % 0: Pencil beam, 1: Isotropically emitting line or point source, 2: Infinite plane wave, 3: Laguerre-Gaussian LG01 beam, 4: Radial-factorizable beam (e.g., a Gaussian beam), 5: X/Y factorizable beam (e.g., a rectangular LED emitter)
model1.MC.lightSource.xFocus       = 0; % [cm] x position of focus
model1.MC.lightSource.yFocus       = -0.01; % [cm] y position of focus
model1.MC.lightSource.zFocus       = model1.G.Lz/2; % [cm] z position of focus

model1.MC.lightSource.theta        = 0; % [rad] Polar angle of beam center axis
model1.MC.lightSource.phi          = 0; % [rad] Azimuthal angle of beam center axis

% Set the power of model1
model1.MC.P                        = 1;

% Copy the model to another object. We now have two identical, but
% separate, model objects
model2 = model1;

% Change model2 to have a different yFocus and a different power than
% model1, but otherwise they are still the same
model2.MC.lightSource.yFocus       = 0.01;
model2.MC.P                        = 0.2;

% Run MC on the two objects individually
model1 = runMonteCarlo(model1);
model2 = runMonteCarlo(model2);

% Combine the two models. combineModels will add up the 
combinedModel = combineModels([model1, model2],'MC'); % First argument is an array of model objects. Second argument is either 'MC' or 'FMC'.


plot(combinedModel,'MC');

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
  mediaProperties(j).mus   = 100; % [cm^-1]
  mediaProperties(j).g     = 0.9;
end