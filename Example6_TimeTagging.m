%% Description
% This example is concerned with time-tagging of the photons. The
% time-of-flight is registered and binned when(if) the photon hits the
% detector, which in this demonstration is spatially-resolved in addition
% to time-resolved. Therefore, the output "image" is 3D, with two spatial
% dimensions and one time-of-flight dimension. The geometry consists of big
% scattering voxels placed diagonally along the xyz direction and along the
% yz direction, illuminated with an infinite plane wave. The xyz-diagonally
% placed voxels are all in the focal plane of the detection lens, so they
% all appear sharp in the time-resolved image, while the yz-diagonally
% placed voxels are not all in the focal plane and some of them are
% therefore blurred out in the image. Scattering from voxels at larger z
% depths are seen to arrive at later times.

%% Geometry definition
model = MCmatlab.model;

model.G.nx                = 20; % Number of bins in the x direction
model.G.ny                = 20; % Number of bins in the y direction
model.G.nz                = 20; % Number of bins in the z direction
model.G.Lx                = .1; % [cm] x size of simulation cuboid
model.G.Ly                = .1; % [cm] y size of simulation cuboid
model.G.Lz                = .1; % [cm] z size of simulation cuboid

model.G.mediaPropertiesFunc = @mediaPropertiesFunc; % Media properties defined as a function at the end of this file
model.G.geomFunc          = @geometryDefinition_TimeTaggingExample; % Function to use for defining the distribution of media in the cuboid. Defined at the end of this m file.

plot(model,'G');

%% Monte Carlo simulation
model.MC.useAllCPUs               = true; % If false, MCmatlab will leave one processor unused. Useful for doing other work on the PC while simulations are running.
model.MC.simulationTimeRequested  = .5; % [min] Time duration of the simulation

model.MC.matchedInterfaces        = true; % Assumes all refractive indices are the same
model.MC.boundaryType             = 1; % 0: No escaping boundaries, 1: All cuboid boundaries are escaping, 2: Top cuboid boundary only is escaping
model.MC.wavelength               = 532; % [nm] Excitation wavelength, used for determination of optical properties for excitation light

model.MC.beam.beamType            = 2; % 0: Pencil beam, 1: Isotropically emitting point source, 2: Infinite plane wave, 3: Laguerre-Gaussian LG01 beam, 4: Radial-factorizable beam (e.g., a Gaussian beam), 5: X/Y factorizable beam (e.g., a rectangular LED emitter)
model.MC.beam.xFocus              = 0; % [cm] x position of focus
model.MC.beam.yFocus              = 0; % [cm] y position of focus
model.MC.beam.zFocus              = model.G.Lz/2; % [cm] z position of focus
model.MC.beam.theta               = 0; % [rad] Polar angle of beam center axis
model.MC.beam.phi                 = 0; % [rad] Azimuthal angle of beam center axis

model.MC.useLightCollector = true;
model.MC.LC.x         = 0; % [cm] x position of either the center of the objective lens focal plane or the fiber tip
model.MC.LC.y         = 0; % [cm] y position
model.MC.LC.z         = model.G.Lz/2; % [cm] z position

model.MC.LC.theta     = atan(1/sqrt(2)); % [rad] Polar angle of direction the light collector is facing
model.MC.LC.phi       = -3*pi/4; % [rad] Azimuthal angle of direction the light collector is facing

model.MC.LC.f         = .2; % [cm] Focal length of the objective lens (if light collector is a fiber, set this to Inf).
model.MC.LC.diam      = .2; % [cm] Diameter of the light collector aperture. For an ideal thin lens, this is 2*f*tan(asin(NA)).
model.MC.LC.fieldSize = .2; % [cm] Field Size of the imaging system (diameter of area in object plane that gets imaged). Only used for finite f.
model.MC.LC.NA        = 0.22; % [-] Fiber NA. Only used for infinite f.

model.MC.LC.res       = 100; % X and Y resolution of light collector in pixels, only used for finite f

model.MC.LC.tStart    = -1.5e-13; % [s] Start of the detection time-of-flight interval
model.MC.LC.tEnd      = 5.5e-12; % [s] End of the detection time-of-flight interval
model.MC.LC.nTimeBins = 100; % Number of bins between tStart and tEnd. If zero, the measurement is not time-resolved.

% Execution, do not modify the next line:
model = runMonteCarlo(model);

plot(model,'MC');

%% Post-processing

%% Geometry function(s)
% A geometry function takes as input X,Y,Z matrices as returned by the
% "ndgrid" MATLAB function as well as any parameters the user may have
% provided in the definition of model.G. It returns the media matrix M,
% containing numerical values indicating the media type (as defined in
% mediaPropertiesFunc) at each voxel location.
function M = geometryDefinition_TimeTaggingExample(X,Y,Z,parameters)
[nx,ny,~] = size(X);
M = ones(size(X)); % Air background
M(1:(nx*(ny+1)+1):end) = 2; % Set xyz diagonal positions to test scatterer
M(1:(nx*(ny+1)):end) = 2; % Set yz diagonal positions to test scatterer
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
mediaProperties(j).mua   = 1e-8;
mediaProperties(j).mus   = 1e-8;
mediaProperties(j).g     = 1;

j=2;
mediaProperties(j).name  = 'test scatterer';
mediaProperties(j).mua   = 0.0000001;
mediaProperties(j).mus   = 100;
mediaProperties(j).g     = 0;
end