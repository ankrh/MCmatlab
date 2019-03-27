addpath([fileparts(matlab.desktop.editor.getActiveFilename) '/helperfuncs']); % The helperfuncs folder is added to the path for the duration of this MATLAB session

%% Description
% This example shows how to use the calcFdet flag to calculate and plot the
% fluence rate of only that light which ends up on the light collector.
% This is shown for both excitation light and fluorescence light. The
% geometry is the same as in example 4, into which a Gaussian beam is
% injected at x = 0.02 and the light collector is looking at x = -0.02.
% 
% In the fluence rate plot for collected light, you can see how the photons
% all start at the source and end at the light collector.

%% Geometry definition
clear Ginput
Ginput.matchedInterfaces = true; % Assumes all refractive indices are 1
Ginput.boundaryType      = 1; % 0: No escaping boundaries, 1: All cuboid boundaries are escaping, 2: Top cuboid boundary only is escaping

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
MCinput.simulationTime           = .5; % [min] Time duration of the simulation
MCinput.calcF                    = true; % (Default: true) If true, the 3D fluence rate output matrix F will be calculated. Set to false if you have a light collector and you're only interested in the Image output.
MCinput.calcFdet                 = true; % (Default: false) If true, the 3D fluence rate output matrix Fdet will be calculated. Only photons that end up on the light collector are counted in Fdet.

MCinput.Beam.beamType            = 3; % 0: Pencil beam, 1: Isotropically emitting point source, 2: Infinite plane wave, 3: Gaussian focus, Gaussian far field beam, 4: Gaussian focus, top-hat far field beam, 5: Top-hat focus, Gaussian far field beam, 6: Top-hat focus, top-hat far field beam, 7: Laguerre-Gaussian LG01 beam
MCinput.Beam.xFocus              = 0.02; % [cm] x position of focus
MCinput.Beam.yFocus              = 0; % [cm] y position of focus
MCinput.Beam.zFocus              = 0; % [cm] z position of focus
MCinput.Beam.theta               = 0; % [rad] Polar angle of beam center axis
MCinput.Beam.phi                 = 0; % [rad] Azimuthal angle of beam center axis
MCinput.Beam.waist               = 0.005; % [cm] Beam waist 1/e^2 radius
MCinput.Beam.divergence          = 5/180*pi; % [rad] Beam divergence 1/e^2 half-angle of beam (for a diffraction limited Gaussian beam, this is G.wavelength*1e-9/(pi*MCinput.Beam.waist*1e-2))

MCinput.LightCollector.x         = -0.02; % [cm] x position of either the center of the objective lens focal plane or the fiber tip
MCinput.LightCollector.y         = 0; % [cm] y position
MCinput.LightCollector.z         = 0; % [cm] z position

MCinput.LightCollector.theta     = 0; % [rad] Polar angle of direction the light collector is facing
MCinput.LightCollector.phi       = pi/2; % [rad] Azimuthal angle of direction the light collector is facing

MCinput.LightCollector.f         = .1; % [cm] Focal length of the objective lens (if light collector is a fiber, set this to Inf).
MCinput.LightCollector.diam      = .1; % [cm] Diameter of the light collector aperture. For an ideal thin lens, this is 2*f*tan(asin(NA)).
MCinput.LightCollector.FieldSize = .04; % [cm] Field Size of the imaging system (diameter of area in object plane that gets imaged). Only used for finite f.
MCinput.LightCollector.NA        = 0.22; % [-] Fiber NA. Only used for infinite f.

MCinput.LightCollector.res       = 50; % X and Y resolution of light collector in pixels, only used for finite f

% Execution, do not modify the next three lines:
MCinput.G = Goutput;
MCoutput = runMonteCarlo(MCinput);
plotMCmatlab(MCinput,MCoutput);

%% Fluorescence Monte Carlo
clear FMCinput
FMCinput.simulationTime           = .1; % [min] Time duration of the simulation
FMCinput.calcF                    = true; % (Default: true) If true, the 3D fluence rate output matrix F will be calculated. Set to false if you have a light collector and you're only interested in the Image output.
FMCinput.calcFdet                 = true; % (Default: false) If true, the 3D fluence rate output matrix Fdet will be calculated. Only photons that end up on the light collector are counted in Fdet.

FMCinput.LightCollector.x         = -0.02; % [cm] x position of either the center of the objective lens focal plane or the fiber tip
FMCinput.LightCollector.y         = 0; % [cm] y position
FMCinput.LightCollector.z         = 0; % [cm] z position

FMCinput.LightCollector.theta     = 0; % [rad] Polar angle of direction the light collector is facing
FMCinput.LightCollector.phi       = pi/2; % [rad] Azimuthal angle of direction the light collector is facing

FMCinput.LightCollector.f         = .1; % [cm] Focal length of the objective lens (if light collector is a fiber, set this to Inf).
FMCinput.LightCollector.diam      = .1; % [cm] Diameter of the light collector aperture. For an ideal thin lens, this is 2*f*tan(asin(NA)).
FMCinput.LightCollector.FieldSize = .04; % [cm] Field Size of the imaging system (diameter of area in object plane that gets imaged). Only used for finite f.
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
