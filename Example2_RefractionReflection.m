addpath([fileparts(matlab.desktop.editor.getActiveFilename) '/helperfuncs']); % The helperfuncs folder is added to the path for the duration of this MATLAB session

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

%% Geometry definition
clear Ginput
Ginput.matchedInterfaces = false; % Uses the refractive indices as defined in getMediaProperties
Ginput.boundaryType      = 1; % 0: No escaping boundaries, 1: All cuboid boundaries are escaping, 2: Top cuboid boundary only is escaping

Ginput.wavelength        = 532; % [nm] Excitation wavelength, used for determination of optical properties for excitation light

Ginput.nx                = 100; % Number of bins in the x direction
Ginput.ny                = 100; % Number of bins in the y direction
Ginput.nz                = 100; % Number of bins in the z direction
Ginput.Lx                = .1; % [cm] x size of simulation cuboid
Ginput.Ly                = .1; % [cm] y size of simulation cuboid
Ginput.Lz                = .1; % [cm] z size of simulation cuboid

Ginput.getCustomMediaProperties = @customMediaProperties; % Custom media properties defined as a function at the bottom of this file

Ginput.GeomFunc          = @GeometryDefinition_RefractionReflectionExample; % Function to use for defining the distribution of media in the cuboid. Defined at the end of this m file.

% Execution, do not modify the next two lines:
Goutput = defineGeometry(Ginput);
plotMCmatlabGeom(Goutput);

%% Monte Carlo simulation
clear MCinput
MCinput.useAllCPUs               = true; % If false, MCmatlab will leave one processor unused. Useful for doing other work on the PC while simulations are running.
MCinput.simulationTime           = .1; % [min] Time duration of the simulation
MCinput.nExamplePaths            = 100; % (Default: 0) This number of photons will have their paths stored and shown after completion, for illustrative purposes

MCinput.Beam.beamType            = 1; % 0: Pencil beam, 1: Isotropically emitting point source, 2: Infinite plane wave, 3: Gaussian focus, Gaussian far field beam, 4: Gaussian focus, top-hat far field beam, 5: Top-hat focus, Gaussian far field beam, 6: Top-hat focus, top-hat far field beam, 7: Laguerre-Gaussian LG01 beam
MCinput.Beam.xFocus              = 0; % [cm] x position of focus
MCinput.Beam.yFocus              = 0; % [cm] y position of focus
MCinput.Beam.zFocus              = 0.04; % [cm] z position of focus
MCinput.Beam.theta               = 0; % [rad] Polar angle of beam center axis
MCinput.Beam.phi                 = 0; % [rad] Azimuthal angle of beam center axis
MCinput.Beam.waist               = 0.005; % [cm] Beam waist 1/e^2 radius
MCinput.Beam.divergence          = 5/180*pi; % [rad] Beam divergence 1/e^2 half-angle of beam (for a diffraction limited Gaussian beam, this is G.wavelength*1e-9/(pi*MCinput.Beam.waist*1e-2))

% Execution, do not modify the next three lines:
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
function M = GeometryDefinition_RefractionReflectionExample(X,Y,Z,parameters)
M = ones(size(X)); % Air background
M(Z>0.03) = 2; % Water
M(Z>0.09) = 3; % Reflector
end

%% Custom media properties
% a function returning struct containing custom media properties for the model.
% For more details on how to define them, check mediaPropertiesLibrary.m
function mediaProperties = customMediaProperties(wavelength,parameters)
j=3;
mediaProperties(j).name  = 'reflector';
mediaProperties(j).mua   = 1;
mediaProperties(j).mus   = 1;
mediaProperties(j).g     = 0;
mediaProperties(j).n     = Inf;
end