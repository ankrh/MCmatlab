addpath([fileparts(matlab.desktop.editor.getActiveFilename) '/helperfuncs']); % The helperfuncs folder is added to the path for the duration of this MATLAB session
fprintf('\n');

%% Decription
% This example showcases the option of letting the thermal media properties
% depend on temperature. A medium is defined that has a volumetric heat
% capacity (VHC) that increases drastically when the temperature is around
% 40 deg C. Such a medium is not particularly realistic, but serves to
% showcase the feature.
% 
% As the medium heats up, the regions that approach 40 deg C will not heat
% much further due to the very high VHC, despite still absorbing the same
% power as before. In the plot, the yellow color (40 deg C) can be seen to
% extend further and further down into the medium as time goes, but even
% the top layer has not become white (45 deg C) by the end of the
% simulation.
% 
% In models where optical or thermal properties depend on temperature or
% fractional damage, it is very important for the user to test whether
% mediaProperties(j).nBins and model.HS.nUpdates are set high enough and
% model.HS.mediaPropRecalcPeriod low enough for a suitably converged
% result.

%% Geometry definition
model = initializeMCmatlabModel();

model.G.nx                = 100; % Number of bins in the x direction
model.G.ny                = 100; % Number of bins in the y direction
model.G.nz                = 100; % Number of bins in the z direction
model.G.Lx                = .1; % [cm] x size of simulation cuboid
model.G.Ly                = .1; % [cm] y size of simulation cuboid
model.G.Lz                = .1; % [cm] z size of simulation cuboid

model.G.mediaPropertiesFunc = @mediaPropertiesFunc; % Media properties defined as a function at the end of this file
model.G.geomFunc          = @geometryDefinition_TdependentVHC; % Function to use for defining the distribution of media in the cuboid. Defined at the end of this m file.

% Execution, do not modify the next line:
model = defineGeometry(model);

plotMCmatlabGeom(model);

%% Monte Carlo simulation
model = clearMCmatlabModel(model,'MC'); % Only necessary if you want to run this section repeatedly, re-using previous G data

model.MC.simulationTime           = .1; % [min] Time duration of the simulation

model.MC.matchedInterfaces        = true; % Assumes all refractive indices are 1
model.MC.boundaryType             = 1; % 0: No escaping boundaries, 1: All cuboid boundaries are escaping, 2: Top cuboid boundary only is escaping
model.MC.wavelength               = 532; % [nm] Excitation wavelength, used for determination of optical properties for excitation light

model.MC.beam.beamType            = 4; % 0: Pencil beam, 1: Isotropically emitting point source, 2: Infinite plane wave, 3: Laguerre-Gaussian LG01 beam, 4: Radial-factorizable beam (e.g., a Gaussian beam), 5: X/Y factorizable beam (e.g., a rectangular LED emitter)
model.MC.beam.NF.radialDistr      = 0; % Radial near field distribution - 0: Top-hat, 1: Gaussian, Array: Custom. Doesn't need to be normalized.
model.MC.beam.NF.radialWidth      = .03; % [cm] Radial near field 1/e^2 radius if top-hat or Gaussian or half-width of the full distribution if custom
model.MC.beam.FF.radialDistr      = 0; % Radial far field distribution - 0: Top-hat, 1: Gaussian, 2: Cosine (Lambertian), Array: Custom. Doesn't need to be normalized.
model.MC.beam.FF.radialWidth      = 0; % [rad] Radial far field 1/e^2 half-angle if top-hat or Gaussian or half-angle of the full distribution if custom. For a diffraction limited Gaussian beam, this should be set to model.MC.wavelength*1e-9/(pi*model.MC.beam.NF.radialWidth*1e-2))
model.MC.beam.xFocus              = 0; % [cm] x position of focus
model.MC.beam.yFocus              = 0; % [cm] y position of focus
model.MC.beam.zFocus              = 0; % [cm] z position of focus
model.MC.beam.theta               = 0; % [rad] Polar angle of beam center axis
model.MC.beam.phi                 = 0; % [rad] Azimuthal angle of beam center axis

% Execution, do not modify the next line:
model = runMonteCarlo(model);

plotMCmatlab(model);

%% Heat simulation
model = clearMCmatlabModel(model,'HS'); % Only necessary if you want to run this section repeatedly, re-using previous G, MC and/or FMC data

model.MC.P                   = 2; % [W] Incident pulse peak power (in case of infinite plane waves, only the power incident upon the cuboid's top surface)

model.HS.largeTimeSteps      = true; % (Default: false) If true, calculations will be faster, but some voxel temperatures may be slightly less precise. Test for yourself whether this precision is acceptable for your application.

model.HS.heatBoundaryType    = 0; % 0: Insulating boundaries, 1: Constant-temperature boundaries (heat-sinked)
model.HS.durationOn          = 0.02; % [s] Pulse on-duration
model.HS.durationOff         = 0.000; % [s] Pulse off-duration
model.HS.durationEnd         = 0.00; % [s] Non-illuminated relaxation time to add to the end of the simulation to let temperature diffuse after the pulse train
model.HS.Tinitial            = 20; % [deg C] Initial temperature

model.HS.plotTempLimits      = [20 45]; % [deg C] Expected range of temperatures, used only for setting the color scale in the plot
model.HS.nUpdates            = 100; % Number of times data is extracted for plots during each pulse. A minimum of 1 update is performed in each phase (2 for each pulse consisting of an illumination phase and a diffusion phase)
model.HS.mediaPropRecalcPeriod = 1; % Every N updates, the media properties will be recalculated (including, if needed, re-running MC and FMC steps)
model.HS.slicePositions      = [.5 0.6 1]; % Relative slice positions [x y z] for the 3D plots on a scale from 0 to 1
model.HS.tempSensorPositions = [0 0 0.005
                                0 0 0.015
                                0 0 0.035
                                0 0 0.065]; % Each row is a temperature sensor's absolute [x y z] coordinates. Leave the matrix empty ([]) to disable temperature sensors.

% Execution, do not modify the next line:
model = simulateHeatDistribution(model);

plotMCmatlabHeat(model);

%% Post-processing

%% Geometry function(s)
% A geometry function takes as input X,Y,Z matrices as returned by the
% "ndgrid" MATLAB function as well as any parameters the user may have
% provided in the definition of Ginput. It returns the media matrix M,
% containing numerical values indicating the media type (as defined in
% mediaPropertiesFunc) at each voxel location.
function M = geometryDefinition_TdependentVHC(X,Y,Z,parameters)
M = ones(size(X)); % fill background with water
M(Z > 0.01) = 2; % Temperature dependent VHC material
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
mediaProperties(j).name  = 'water';
mediaProperties(j).mua   = 0.00036;
mediaProperties(j).mus   = 10;
mediaProperties(j).g     = 1.0;
mediaProperties(j).VHC   = 4.19;
mediaProperties(j).TC    = 5.8e-3;

j=2;
mediaProperties(j).name  = 'temp. dep. VHC material';
mediaProperties(j).mua = 30;
mediaProperties(j).mus = 50;
mediaProperties(j).g   = 0.8;
mediaProperties(j).VHC = '3.5 + 500./(1+exp(40-T))';
mediaProperties(j).TC  = 0.52e-2;
mediaProperties(j).nBins = 30;
end