addpath([fileparts(matlab.desktop.editor.getActiveFilename) '/helperfuncs']); % The helperfuncs folder is added to the path for the duration of this MATLAB session
fprintf('\n');

%% Description
% This example shows how to use the calcNFRdet flag to calculate and plot
% the fluence rate of only that light which ends up on the light collector.
% This is shown for both excitation light and fluorescence light. The
% geometry is again almost the same as in example 5, into which a Gaussian
% beam is injected at x = 0.02 and the light collector is looking at x =
% -0.02.
% 
% In the fluence rate plot for collected light, you can see how the photons
% all start at the source and end at the light collector.
% 
% This example also shows two other features: (1) That the Monte Carlo
% simulation can be set to launch a set number of photons rather than run
% for a set time using the MC.nPhotonsRequested and FMC.nPhotonsRequested
% fields, and (2) that by setting the boundaryType flags to 2, we can allow
% photons to travel outside the cuboid in the x, y, and +z directions as
% seen in the photon paths illustrations, although absorption and fluence
% rate data is still only tracked within the main cuboid.

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

model.MC.nPhotonsRequested        = 5e6; % # of photons to launch

model.MC.calcNFR                  = true; % (Default: true) If true, the 3D fluence rate output matrix NFR will be calculated. Set to false if you have a light collector and you're only interested in the image output.
model.MC.calcNFRdet               = true; % (Default: false) If true, the 3D fluence rate output matrix NFRdet will be calculated. Only photons that end up on the light collector are counted in NFRdet.
model.MC.nExamplePaths            = 200;

model.MC.matchedInterfaces        = true; % Assumes all refractive indices are 1
model.MC.boundaryType             = 2; % 0: No escaping boundaries, 1: All cuboid boundaries are escaping, 2: Top cuboid boundary only is escaping
model.MC.wavelength               = 450; % [nm] Excitation wavelength, used for determination of optical properties for excitation light

model.MC.beam.beamType            = 4; % 0: Pencil beam, 1: Isotropically emitting point source, 2: Infinite plane wave, 3: Laguerre-Gaussian LG01 beam, 4: Radial-factorizable beam (e.g., a Gaussian beam), 5: X/Y factorizable beam (e.g., a rectangular LED emitter)
model.MC.beam.NF.radialDistr      = 1; % Radial near field distribution - 0: Top-hat, 1: Gaussian, Array: Custom. Doesn't need to be normalized.
model.MC.beam.NF.radialWidth      = .005; % [cm] Radial near field 1/e^2 radius if top-hat or Gaussian or half-width of the full distribution if custom
model.MC.beam.FF.radialDistr      = 1; % Radial far field distribution - 0: Top-hat, 1: Gaussian, 2: Cosine (Lambertian), Array: Custom. Doesn't need to be normalized.
model.MC.beam.FF.radialWidth      = 5/180*pi; % [rad] Radial far field 1/e^2 half-angle if top-hat or Gaussian or half-angle of the full distribution if custom. For a diffraction limited Gaussian beam, this should be set to model.MC.wavelength*1e-9/(pi*model.MC.beam.NF.radialWidth*1e-2))
model.MC.beam.xFocus              = 0.02; % [cm] x position of focus
model.MC.beam.yFocus              = 0; % [cm] y position of focus
model.MC.beam.zFocus              = 0; % [cm] z position of focus
model.MC.beam.theta               = 0; % [rad] Polar angle of beam center axis
model.MC.beam.phi                 = 0; % [rad] Azimuthal angle of beam center axis

model.MC.useLightCollector        = true;
model.MC.LC.x                     = -0.02; % [cm] x position of either the center of the objective lens focal plane or the fiber tip
model.MC.LC.y                     = 0; % [cm] y position
model.MC.LC.z                     = 0; % [cm] z position

model.MC.LC.theta                 = 0; % [rad] Polar angle of direction the light collector is facing
model.MC.LC.phi                   = pi/2; % [rad] Azimuthal angle of direction the light collector is facing

model.MC.LC.f                     = .1; % [cm] Focal length of the objective lens (if light collector is a fiber, set this to Inf).
model.MC.LC.diam                  = .1; % [cm] Diameter of the light collector aperture. For an ideal thin lens, this is 2*f*tan(asin(NA)).
model.MC.LC.fieldSize             = .04; % [cm] Field Size of the imaging system (diameter of area in object plane that gets imaged). Only used for finite f.
model.MC.LC.NA                    = 0.22; % [-] Fiber NA. Only used for infinite f.

model.MC.LC.res                   = 50; % X and Y resolution of light collector in pixels, only used for finite f

% Execution, do not modify the next line:
model = runMonteCarlo(model);

plotMCmatlab(model);

%% Fluorescence Monte Carlo
model = clearMCmatlabModel(model,'FMC'); % Only necessary if you want to run this section repeatedly, re-using previous G and MC data

model.FMC.nPhotonsRequested       = 5e6; % # of photons to launch

model.FMC.calcNFR                 = true; % (Default: true) If true, the 3D fluence rate output matrix NFR will be calculated. Set to false if you have a light collector and you're only interested in the image output.
model.FMC.calcNFRdet              = true; % (Default: false) If true, the 3D fluence rate output matrix NFRdet will be calculated. Only photons that end up on the light collector are counted in NFRdet.
model.FMC.nExamplePaths           = 200;

model.FMC.matchedInterfaces       = true; % Assumes all refractive indices are 1
model.FMC.boundaryType            = 2; % 0: No escaping boundaries, 1: All cuboid boundaries are escaping, 2: Top cuboid boundary only is escaping
model.FMC.wavelength              = 550; % [nm] Fluorescence wavelength, used for determination of optical properties for fluorescence light

model.FMC.useLightCollector       = true;
model.FMC.LC.x                    = -0.02; % [cm] x position of either the center of the objective lens focal plane or the fiber tip
model.FMC.LC.y                    = 0; % [cm] y position
model.FMC.LC.z                    = 0; % [cm] z position

model.FMC.LC.theta                = 0; % [rad] Polar angle of direction the light collector is facing
model.FMC.LC.phi                  = pi/2; % [rad] Azimuthal angle of direction the light collector is facing

model.FMC.LC.f                    = .1; % [cm] Focal length of the objective lens (if light collector is a fiber, set this to Inf).
model.FMC.LC.diam                 = .1; % [cm] Diameter of the light collector aperture. For an ideal thin lens, this is 2*f*tan(asin(NA)).
model.FMC.LC.fieldSize            = .04; % [cm] Field Size of the imaging system (diameter of area in object plane that gets imaged). Only used for finite f.
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
% fractional heat damage FD can be specified as in examples 12-15.
function mediaProperties = mediaPropertiesFunc(wavelength,parameters)
j=1;
mediaProperties(j).name  = 'fluorescence absorber';
if(wavelength<500)
  mediaProperties(j).mua = 10;
  mediaProperties(j).mus = 100;
  mediaProperties(j).g   = 0.9;
else
  mediaProperties(j).mua = 100;
  mediaProperties(j).mus = 100;
  mediaProperties(j).g   = 0.9;
end
mediaProperties(j).n   = 1.3;

j=2;
mediaProperties(j).name  = 'fluorescer';
if(wavelength<500)
  mediaProperties(j).mua = 100;
  mediaProperties(j).mus = 100;
  mediaProperties(j).g   = 0.9;
  
  % Only one of PY and QY may be defined:
  mediaProperties(j).PY   = 0.5; % Fluorescence power yield (ratio of power emitted to power absorbed)
%   mediaProperties(j).QY   = 0.6; % Fluorescence quantum yield (ratio of photons emitted to photons absorbed)
else
  mediaProperties(j).mua = 10;
  mediaProperties(j).mus = 100;
  mediaProperties(j).g   = 0.9;
end
mediaProperties(j).n   = 1.3;

end
