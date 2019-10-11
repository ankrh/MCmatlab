addpath([fileparts(matlab.desktop.editor.getActiveFilename) '/helperfuncs']); % The helperfuncs folder is added to the path for the duration of this MATLAB session

%% Description
% 

%% Geometry definition
clear Ginput
Ginput.nx                = 100; % Number of bins in the x direction
Ginput.ny                = 100; % Number of bins in the y direction
Ginput.nz                = 100; % Number of bins in the z direction
Ginput.Lx                = .1; % [cm] x size of simulation cuboid
Ginput.Ly                = .1; % [cm] y size of simulation cuboid
Ginput.Lz                = .1; % [cm] z size of simulation cuboid

Ginput.mediaPropertiesFunc = @mediaPropertiesFunc; % Media properties defined as a function at the end of this file
Ginput.GeomFunc          = @GeometryDefinition_SaturableAbsorber; % Function to use for defining the distribution of media in the cuboid. Defined at the end of this m file.

% Execution, do not modify the next line:
Goutput = defineGeometry(Ginput);

plotMCmatlabGeom(Goutput);

%% Monte Carlo simulation
clear MCinput
MCinput.simulationTime           = .1; % [min] Time duration of the simulation

MCinput.matchedInterfaces        = true; % Assumes all refractive indices are 1
MCinput.boundaryType             = 1; % 0: No escaping boundaries, 1: All cuboid boundaries are escaping, 2: Top cuboid boundary only is escaping
MCinput.wavelength               = 532; % [nm] Excitation wavelength, used for determination of optical properties for excitation light

MCinput.Beam.beamType            = 3; % 0: Pencil beam, 1: Isotropically emitting point source, 2: Infinite plane wave, 3: Gaussian focus, Gaussian far field beam, 4: Gaussian focus, top-hat far field beam, 5: Top-hat focus, Gaussian far field beam, 6: Top-hat focus, top-hat far field beam, 7: Laguerre-Gaussian LG01 beam
MCinput.Beam.xFocus              = 0; % [cm] x position of focus
MCinput.Beam.yFocus              = 0; % [cm] y position of focus
MCinput.Beam.zFocus              = Ginput.Lz/2; % [cm] z position of focus
MCinput.Beam.theta               = 0; % [rad] Polar angle of beam center axis
MCinput.Beam.phi                 = 0; % [rad] Azimuthal angle of beam center axis
MCinput.Beam.waist               = 0.02; % [cm] Beam waist 1/e^2 radius
MCinput.Beam.divergence          = 0; % [rad] Beam divergence 1/e^2 half-angle of beam (for a diffraction limited Gaussian beam, this is G.wavelength*1e-9/(pi*MCinput.Beam.waist*1e-2))

MCinput.Iinitial = rand(Ginput.nx,Ginput.ny,Ginput.nz);

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
function M = GeometryDefinition_SaturableAbsorber(X,Y,Z,parameters)
absorberdepth = 0.03;
M = ones(size(X)); % Air
M(Z > absorberdepth) = 2; % Saturable absorber
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
mediaProperties(j).name  = 'saturable absorber';
mediaProperties(j).mua   = '1./(1+I)';
mediaProperties(j).mus   = 1;
mediaProperties(j).g     = 0.9;
end