addpath([fileparts(matlab.desktop.editor.getActiveFilename) '/helperfuncs']); % The helperfuncs folder is added to the path for the duration of this MATLAB session
fprintf('\n');

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
model = initializeMCmatlabModel();

model.G.nx                = 100; % Number of bins in the x direction
model.G.ny                = 100; % Number of bins in the y direction
model.G.nz                = 100; % Number of bins in the z direction
model.G.Lx                = .1; % [cm] x size of simulation cuboid
model.G.Ly                = .1; % [cm] y size of simulation cuboid
model.G.Lz                = .1; % [cm] z size of simulation cuboid

model.G.mediaPropertiesFunc = @mediaPropertiesFunc; % Media properties defined as a function at the end of this file
model.G.geomFunc          = @geometryDefinition_FluorescingCylinder; % Function to use for defining the distribution of media in the cuboid. Defined at the end of this m file.

% Execution, do not modify the next line:
model = defineGeometry(model);

plotMCmatlabGeom(model);

%% Monte Carlo simulation
model = clearMCmatlabModel(model,'MC'); % Only necessary if you want to run this section repeatedly, re-using previous G data

model.MC.useAllCPUs               = true; % If false, MCmatlab will leave one processor unused. Useful for doing other work on the PC while simulations are running.
model.MC.simulationTime           = .1; % [min] Time duration of the simulation
model.MC.nExamplePaths            = 100; % (Default: 0) This number of photons will have their paths stored and shown after completion, for illustrative purposes

model.MC.matchedInterfaces        = true; % Assumes all refractive indices are 1
model.MC.boundaryType             = 1; % 0: No escaping boundaries, 1: All cuboid boundaries are escaping, 2: Top cuboid boundary only is escaping
model.MC.wavelength               = 450; % [nm] Excitation wavelength, used for determination of optical properties for excitation light

model.MC.beam.beamType            = 2; % 0: Pencil beam, 1: Isotropically emitting point source, 2: Infinite plane wave, 3: Gaussian focus, Gaussian far field beam, 4: Gaussian focus, top-hat far field beam, 5: Top-hat focus, Gaussian far field beam, 6: Top-hat focus, top-hat far field beam, 7: Laguerre-Gaussian LG01 beam
model.MC.beam.xFocus              = 0; % [cm] x position of focus
model.MC.beam.yFocus              = 0; % [cm] y position of focus
model.MC.beam.zFocus              = model.G.Lz/2; % [cm] z position of focus
model.MC.beam.theta               = 0; % [rad] Polar angle of beam center axis
model.MC.beam.phi                 = 0; % [rad] Azimuthal angle of beam center axis
model.MC.beam.waist               = 0.005; % [cm] Beam waist 1/e^2 radius
model.MC.beam.divergence          = 5/180*pi; % [rad] Beam divergence 1/e^2 half-angle of beam (for a diffraction limited Gaussian beam, this is G.wavelength*1e-9/(pi*model.MC.beam.waist*1e-2))

model.MC.useLightCollector        = true;

model.MC.LC.x                     = 0; % [cm] x position of either the center of the objective lens focal plane or the fiber tip
model.MC.LC.y                     = 0; % [cm] y position
model.MC.LC.z                     = 0.03; % [cm] z position

model.MC.LC.theta                 = 0; % [rad] Polar angle of direction the light collector is facing
model.MC.LC.phi                   = pi/2; % [rad] Azimuthal angle of direction the light collector is facing

model.MC.LC.f                     = .2; % [cm] Focal length of the objective lens (if light collector is a fiber, set this to Inf).
model.MC.LC.diam                  = .1; % [cm] Diameter of the light collector aperture. For an ideal thin lens, this is 2*f*tan(asin(NA)).
model.MC.LC.fieldSize             = .1; % [cm] Field Size of the imaging system (diameter of area in object plane that gets imaged). Only used for finite f.
model.MC.LC.NA                    = 0.22; % [-] Fiber NA. Only used for infinite f.

model.MC.LC.res                   = 50; % X and Y resolution of light collector in pixels, only used for finite f

% Execution, do not modify the next line:
model = runMonteCarlo(model);

plotMCmatlab(model);

%% Fluorescence Monte Carlo
model = clearMCmatlabModel(model,'FMC'); % Only necessary if you want to run this section repeatedly, re-using previous G and MC data

model.FMC.useAllCPUs              = true; % If false, MCmatlab will leave one processor unused. Useful for doing other work on the PC while simulations are running.
model.FMC.simulationTime          = .1; % [min] Time duration of the simulation
model.FMC.nExamplePaths           = 100; % (Default: 0) This number of photons will have their paths stored and shown after completion, for illustrative purposes

model.FMC.matchedInterfaces       = true; % Assumes all refractive indices are 1
model.FMC.boundaryType            = 1; % 0: No escaping boundaries, 1: All cuboid boundaries are escaping, 2: Top cuboid boundary only is escaping
model.FMC.wavelength              = 550; % [nm] Fluorescence wavelength, used for determination of optical properties for fluorescence light

model.FMC.useLightCollector       = true;

model.FMC.LC.x                    = 0; % [cm] x position of either the center of the objective lens focal plane or the fiber tip
model.FMC.LC.y                    = 0; % [cm] y position
model.FMC.LC.z                    = 0.03; % [cm] z position

model.FMC.LC.theta                = 0; % [rad] Polar angle of direction the light collector is facing
model.FMC.LC.phi                  = pi/2; % [rad] Azimuthal angle of direction the light collector is facing

model.FMC.LC.f                    = .2; % [cm] Focal length of the objective lens (if light collector is a fiber, set this to Inf).
model.FMC.LC.diam                 = .1; % [cm] Diameter of the light collector aperture. For an ideal thin lens, this is 2*f*tan(asin(NA)).
model.FMC.LC.fieldSize            = .1; % [cm] Field Size of the imaging system (diameter of area in object plane that gets imaged). Only used for finite f.
model.FMC.LC.NA                   = 0.22; % [-] Fiber NA. Only used for infinite f.

model.FMC.LC.res                  = 50; % X and Y resolution of light collector in pixels, only used for finite f

% Execution, do not modify the next line:
model = runMonteCarlo(model,'fluorescence');

plotMCmatlab(model,'fluorescence');

%% Post-processing

%% Geometry function(s)
% A geometry function takes as input X,Y,Z matrices as returned by the
% "ndgrid" MATLAB function as well as any parameters the user may have
% provided in the definition of Ginput. It returns the media matrix M,
% containing numerical values indicating the media type (as defined in
% mediaPropertiesFunc) at each voxel location.
function M = geometryDefinition_FluorescingCylinder(X,Y,Z,parameters)
cylinderradius  = 0.0100;
M = ones(size(X)); % fill background with fluorescence absorber
M(Y.^2 + (Z - 3*cylinderradius).^2 < cylinderradius^2) = 2; % fluorescer
end

%% Media Properties function
% The media properties function defines all the optical and thermal
% properties of the media involved by constructing and returning a
% "mediaProperties" struct with various fields. As its input, the function
% takes the wavelength as well as any other parameters you might specify
% above in the model file, for example parameters that you might loop over
% in a for loop. Dependence on excitation fluence rate FR, temperature T or
% fractional heat damage FD can be specified as in examples 11-14.
function mediaProperties = mediaPropertiesFunc(wavelength,parameters)
j=1;
mediaProperties(j).name  = 'test fluorescence absorber';
if(wavelength<500)
  mediaProperties(j).mua = 1;
  mediaProperties(j).mus = 100;
  mediaProperties(j).g   = 0.9;
else
  mediaProperties(j).mua = 100;
  mediaProperties(j).mus = 100;
  mediaProperties(j).g   = 0.9;
end

j=2;
mediaProperties(j).name  = 'test fluorescer';
if(wavelength<500)
  mediaProperties(j).mua = 100;
  mediaProperties(j).mus = 100;
  mediaProperties(j).g   = 0.9;

  mediaProperties(j).Y   = 0.5;
else
  mediaProperties(j).mua = 1;
  mediaProperties(j).mus = 100;
  mediaProperties(j).g   = 0.9;
end
end
