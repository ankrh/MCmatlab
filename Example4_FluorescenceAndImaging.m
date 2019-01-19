addpath([fileparts(matlab.desktop.editor.getActiveFilename) '/helperfuncs']); % The helperfuncs folder is added to the path for the duration of this MATLAB session

%% Description
% In this example, simulation of fluorescence (luminescence) is shown. The
% test geometry is a fluorescing cylinder in which excitation light is
% predominantly absorbed embedded in a block of medium in which
% fluorescence light is predominantly absorbed. The geometry is illuminated
% with an infinite plane wave, for which the xFocus, yFocus, zFocus, waist
% and divergence quantities are not used.
%
% This example also shows detection of the light exiting the cuboid,
% separately for excitation light and for fluorescence light. Although most
% of the fluorescence light is absorbed in the medium surrounding the
% cylinder, some of it escapes to the detector, showing a slightly blurred
% image of the cylinder.
% 
% The "nExamplePaths" parameter (see example 2) is used for both the
% excitation and fluorescence simulations, showing paths of both kinds of
% photons.

%% Geometry definition
clear Ginput
Ginput.matchedInterfaces = true; % Assumes all refractive indices are 1
Ginput.boundaryType      = 1; % 0: No boundaries, 1: All cuboid boundaries, 2: Top cuboid boundary only

Ginput.wavelength        = 450; % [nm] Excitation wavelength, used for determination of optical properties for excitation light
Ginput.wavelength_f      = 550; % [nm] Fluorescence wavelength, used for determination of optical properties for fluorescence light

Ginput.nx                = 100; % Number of bins in the x direction
Ginput.ny                = 100; % Number of bins in the y direction
Ginput.nz                = 100; % Number of bins in the z direction
Ginput.Lx                = .1; % [cm] x size of simulation cuboid
Ginput.Ly                = .1; % [cm] y size of simulation cuboid
Ginput.Lz                = .1; % [cm] z size of simulation cuboid

Ginput.GeomFunc          = @GeometryDefinition_FluorescingCylinder; % Function to use for defining the distribution of media in the cuboid. Defined at the end of this m file.

% Execution, do not modify the next two lines:
Goutput = defineGeometry(Ginput);
plotMCmatlabGeom(Goutput);

%% Monte Carlo simulation
clear MCinput
MCinput.useAllCPUs               = true; % If false, MCmatlab will leave one processor unused. Useful for doing other work on the PC while simulations are running.
MCinput.simulationTime           = .1; % [min] Time duration of the simulation
MCinput.nExamplePaths            = 100; % (Default: 0) This number of photons will have their paths stored and shown after completion, for illustrative purposes

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
MCinput.LightCollector.z         = 0.03; % [cm] z position

MCinput.LightCollector.theta     = 0; % [rad] Polar angle of direction the light collector is facing
MCinput.LightCollector.phi       = pi/2; % [rad] Azimuthal angle of direction the light collector is facing

MCinput.LightCollector.f         = .2; % [cm] Focal length of the objective lens (if light collector is a fiber, set this to Inf).
MCinput.LightCollector.diam      = .1; % [cm] Diameter of the light collector aperture. For an ideal thin lens, this is 2*f*tan(asin(NA)).
MCinput.LightCollector.FieldSize = .1; % [cm] Field Size of the imaging system (diameter of area in object plane that gets imaged). Only used for finite f.
MCinput.LightCollector.NA        = 0.22; % [-] Fiber NA. Only used for infinite f.

MCinput.LightCollector.res       = 50; % X and Y resolution of light collector in pixels, only used for finite f

% Execution, do not modify the next three lines:
MCinput.G = Goutput;
MCoutput = runMonteCarlo(MCinput);
plotMCmatlab(MCinput,MCoutput);

%% Fluorescence Monte Carlo
clear FMCinput
FMCinput.useAllCPUs               = true; % If false, MCmatlab will leave one processor unused. Useful for doing other work on the PC while simulations are running.
FMCinput.simulationTime           = .1; % [min] Time duration of the simulation
FMCinput.nExamplePaths            = 100; % (Default: 0) This number of photons will have their paths stored and shown after completion, for illustrative purposes

FMCinput.LightCollector.x         = 0; % [cm] x position of either the center of the objective lens focal plane or the fiber tip
FMCinput.LightCollector.y         = 0; % [cm] y position
FMCinput.LightCollector.z         = 0.03; % [cm] z position

FMCinput.LightCollector.theta     = 0; % [rad] Polar angle of direction the light collector is facing
FMCinput.LightCollector.phi       = pi/2; % [rad] Azimuthal angle of direction the light collector is facing

FMCinput.LightCollector.f         = .2; % [cm] Focal length of the objective lens (if light collector is a fiber, set this to Inf).
FMCinput.LightCollector.diam      = .1; % [cm] Diameter of the light collector aperture. For an ideal thin lens, this is 2*f*tan(asin(NA)).
FMCinput.LightCollector.FieldSize = .1; % [cm] Field Size of the imaging system (diameter of area in object plane that gets imaged). Only used for finite f.
FMCinput.LightCollector.NA        = 0.22; % [-] Fiber NA. Only used for infinite f.

FMCinput.LightCollector.res       = 50; % X and Y resolution of light collector in pixels, only used for finite f

% Execution, do not modify the next four lines:
FMCinput.G = Goutput;
FMCinput.MCoutput = MCoutput;
FMCoutput = runMonteCarloFluorescence(FMCinput);
plotMCmatlabFluorescence(FMCinput,FMCoutput);

%% Post-processing

%% Geometry function(s)
% A geometry function takes as input X,Y,Z matrices as returned by the
% "ndgrid" MATLAB function as well as any parameters the user may have
% provided in the definition of Ginput. It returns the media matrix M,
% containing numerical values indicating the media type (as defined in
% getMediaProperties) at each voxel location.
function M = GeometryDefinition_FluorescingCylinder(X,Y,Z,parameters)
cylinderradius  = 0.0100;
M = 17*ones(size(X)); % fill background with fluorescence absorber
M(Y.^2 + (Z - 3*cylinderradius).^2 < cylinderradius^2) = 16; % fluorescer
end
