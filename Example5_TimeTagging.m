addpath([fileparts(matlab.desktop.editor.getActiveFilename) '/helperfuncs']); % The helperfuncs folder is added to the path for the duration of this MATLAB session

%% Description
% This example is concerned with time-tagging of the photons. The time of
% arrival is registered and binned when(if) the photon hits the detector,
% which in this demonstration is spatially-resolved in addition to
% time-resolved. Therefore, the output "image" is 3D, with two spatial
% dimensions and one time dimension. The geometry consists of big
% scattering voxels placed diagonally along the xyz direction and along the
% yz direction, illuminated with an infinite plane wave. The xyz-diagonally
% placed voxels are all in the focal plane of the detection lens, so they
% all appear sharp in the time-resolved image, while the yz-diagonally
% placed voxels are not all in the focal plane and some of them are
% therefore blurred out in the image. Scattering from voxels at larger z
% depths are seen to arrive at later times.

%% Geometry definition
clear Ginput
Ginput.nx                = 20; % Number of bins in the x direction
Ginput.ny                = 20; % Number of bins in the y direction
Ginput.nz                = 20; % Number of bins in the z direction
Ginput.Lx                = .1; % [cm] x size of simulation cuboid
Ginput.Ly                = .1; % [cm] y size of simulation cuboid
Ginput.Lz                = .1; % [cm] z size of simulation cuboid

Ginput.mediaPropertiesFunc = @mediaPropertiesFunc; % Media properties defined as a function at the end of this file
Ginput.GeomFunc          = @GeometryDefinition_TimeTaggingExample; % Function to use for defining the distribution of media in the cuboid. Defined at the end of this m file.

% Execution, do not modify the next line:
Goutput = defineGeometry(Ginput);

plotMCmatlabGeom(Goutput);

%% Monte Carlo simulation
clear MCinput
MCinput.useAllCPUs               = true; % If false, MCmatlab will leave one processor unused. Useful for doing other work on the PC while simulations are running.
MCinput.simulationTime           = .5; % [min] Time duration of the simulation

MCinput.matchedInterfaces        = true; % Assumes all refractive indices are 1
MCinput.boundaryType             = 1; % 0: No escaping boundaries, 1: All cuboid boundaries are escaping, 2: Top cuboid boundary only is escaping
MCinput.wavelength               = 532; % [nm] Excitation wavelength, used for determination of optical properties for excitation light

MCinput.Beam.beamType            = 2; % 0: Pencil beam, 1: Isotropically emitting point source, 2: Infinite plane wave, 3: Gaussian focus, Gaussian far field beam, 4: Gaussian focus, top-hat far field beam, 5: Top-hat focus, Gaussian far field beam, 6: Top-hat focus, top-hat far field beam, 7: Laguerre-Gaussian LG01 beam
MCinput.Beam.xFocus              = 0; % [cm] x position of focus
MCinput.Beam.yFocus              = 0; % [cm] y position of focus
MCinput.Beam.zFocus              = Ginput.Lz/2; % [cm] z position of focus
MCinput.Beam.theta               = 0; % [rad] Polar angle of beam center axis
MCinput.Beam.phi                 = 0; % [rad] Azimuthal angle of beam center axis
MCinput.Beam.waist               = 0.005; % [cm] Beam waist 1/e^2 radius
MCinput.Beam.divergence          = 5/180*pi; % [rad] Beam divergence 1/e^2 half-angle of beam (for a diffraction limited Gaussian beam, this is G.wavelength*1e-9/(pi*MCinput.Beam.waist*1e-2))

MCinput.LightCollector.x         = 0; % [cm] x position of either the center of the objective lens focal plane or the fiber tip
MCinput.LightCollector.y         = 0; % [cm] y position
MCinput.LightCollector.z         = Ginput.Lz/2; % [cm] z position

MCinput.LightCollector.theta     = atan(1/sqrt(2)); % [rad] Polar angle of direction the light collector is facing
MCinput.LightCollector.phi       = -3*pi/4; % [rad] Azimuthal angle of direction the light collector is facing

MCinput.LightCollector.f         = .2; % [cm] Focal length of the objective lens (if light collector is a fiber, set this to Inf).
MCinput.LightCollector.diam      = .2; % [cm] Diameter of the light collector aperture. For an ideal thin lens, this is 2*f*tan(asin(NA)).
MCinput.LightCollector.FieldSize = .2; % [cm] Field Size of the imaging system (diameter of area in object plane that gets imaged). Only used for finite f.
MCinput.LightCollector.NA        = 0.22; % [-] Fiber NA. Only used for infinite f.

MCinput.LightCollector.res       = 100; % X and Y resolution of light collector in pixels, only used for finite f

MCinput.LightCollector.tStart    = -1.5e-13; % [s] Start of the detection time interval
MCinput.LightCollector.tEnd      = 5.5e-12; % [s] End of the detection time interval
MCinput.LightCollector.nTimeBins = 100; % Number of bins between tStart and tEnd. If zero, the measurement is not time-resolved.

% Execution, do not modify the next two lines:
MCinput.G = Goutput;
MCoutput = runMonteCarlo(MCinput);

plotMCmatlab(MCinput,MCoutput);

%% Post-processing

%% Geometry function(s)
% A geometry function takes as input X,Y,Z matrices as returned by the
% "ndgrid" MATLAB function as well as any parameters the user may have
% provided in the definition of Ginput. It returns the media matrix M,
% containing numerical values indicating the media type (as defined in
% getMediaProperties) at each voxel location.
function M = GeometryDefinition_TimeTaggingExample(X,Y,Z,parameters)
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
% in a for loop.
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