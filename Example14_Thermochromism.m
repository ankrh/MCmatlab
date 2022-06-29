%% Decription
% This example demonstrates simulations in which the absorption coefficient
% is dependent on temperature, such as in a thermochromic material. For a
% introduction to how fluence rate or temperature dependent optical or
% thermal properties are binned, see example 12.
% 
% The geometry is a wedge of thermochromic material suspended in water over
% a plane sheet of the same thermochromic material. A focused Gaussian beam
% is incident on the wedge with a slight tilt in the -y direction. The mua
% property of the thermochromic material is defined as a char array using
% standard MATLAB syntax in such a way that mua = 150 for T = 20 deg C but
% drops linearly with temperature until it bottoms out at 1e-8 for T = 40
% deg C.
% 
% The light is on for 200 ms and 50 plotting updates are requested.
% HS.mediaPropRecalcPeriod is set to 2 which means that every 2 updates (8
% ms), the heat simulation is paused while the media properties of the
% incident light are recalculated and the fluence rate is found again by
% re-running the Monte Carlo step. Monte Carlo results are plotted both
% before and after the heat simulations.
% 
% In models where optical or thermal properties depend on temperature or
% fractional damage, it is very important for the user to test whether
% mediaProperties(j).nBins and model.HS.nUpdates are set high enough and
% model.HS.mediaPropRecalcPeriod low enough for a suitably converged
% result.

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
model = MCmatlab.model;

model.G.nx                = 100; % Number of bins in the x direction
model.G.ny                = 100; % Number of bins in the y direction
model.G.nz                = 100; % Number of bins in the z direction
model.G.Lx                = .1; % [cm] x size of simulation cuboid
model.G.Ly                = .1; % [cm] y size of simulation cuboid
model.G.Lz                = .1; % [cm] z size of simulation cuboid

model.G.mediaPropertiesFunc = @mediaPropertiesFunc; % Media properties defined as a function at the end of this file
model.G.geomFunc          = @geometryDefinition; % Function to use for defining the distribution of media in the cuboid. Defined at the end of this m file.

plot(model,'G');

%% Monte Carlo simulation
model.MC.simulationTimeRequested  = .02; % [min] Time duration of the simulation

model.MC.matchedInterfaces        = true; % Assumes all refractive indices are the same
model.MC.boundaryType             = 1; % 0: No escaping boundaries, 1: All cuboid boundaries are escaping, 2: Top cuboid boundary only is escaping, 3: Top and bottom boundaries are escaping, while the side boundaries are cyclic
model.MC.wavelength               = 532; % [nm] Excitation wavelength, used for determination of optical properties for excitation light

model.MC.lightSource.sourceType   = 4; % 0: Pencil beam, 1: Isotropically emitting line or point source, 2: Infinite plane wave, 3: Laguerre-Gaussian LG01 beam, 4: Radial-factorizable beam (e.g., a Gaussian beam), 5: X/Y factorizable beam (e.g., a rectangular LED emitter)
model.MC.lightSource.focalPlaneIntensityDistribution.radialDistr = 1; % Radial focal plane intensity distribution - 0: Top-hat, 1: Gaussian, Array: Custom. Doesn't need to be normalized.
model.MC.lightSource.focalPlaneIntensityDistribution.radialWidth = .005; % [cm] Radial focal plane 1/e^2 radius if top-hat or Gaussian or half-width of the full distribution if custom
model.MC.lightSource.angularIntensityDistribution.radialDistr = 1; % Radial angular intensity distribution - 0: Top-hat, 1: Gaussian, 2: Cosine (Lambertian), Array: Custom. Doesn't need to be normalized.
model.MC.lightSource.angularIntensityDistribution.radialWidth = 10/180*pi; % [rad] Radial angular 1/e^2 half-angle if top-hat or Gaussian or half-angle of the full distribution if custom. For a diffraction limited Gaussian beam, this should be set to model.MC.wavelength*1e-9/(pi*model.MC.lightSource.focalPlaneIntensityDistribution.radialWidth*1e-2))
model.MC.lightSource.xFocus       = 0; % [cm] x position of focus
model.MC.lightSource.yFocus       = -0.017; % [cm] y position of focus
model.MC.lightSource.zFocus       = 0.078; % [cm] z position of focus
model.MC.lightSource.theta        = pi/6; % [rad] Polar angle of beam center axis
model.MC.lightSource.phi          = -pi/2; % [rad] Azimuthal angle of beam center axis

model.HS.Tinitial                 = 20; % [deg C] Initial temperature

% Execution, do not modify the next line:
model = runMonteCarlo(model);

plot(model,'MC');

%% Heat simulation
model.MC.P                   = 2; % [W] Incident pulse peak power (in case of infinite plane waves, only the power incident upon the cuboid's top surface)

model.HS.largeTimeSteps      = true; % (Default: false) If true, calculations will be faster, but some voxel temperatures may be slightly less precise. Test for yourself whether this precision is acceptable for your application.

model.HS.heatBoundaryType    = 0; % 0: Insulating boundaries, 1: Constant-temperature boundaries (heat-sinked)
model.HS.durationOn          = 0.2; % [s] Pulse on-duration
model.HS.durationOff         = 0.00; % [s] Pulse off-duration
model.HS.durationEnd         = 0.00; % [s] Non-illuminated relaxation time to add to the end of the simulation to let temperature diffuse after the pulse train

model.HS.plotTempLimits      = [20 50]; % [deg C] Expected range of temperatures, used only for setting the color scale in the plot
model.HS.nUpdates            = 50; % Number of times data is extracted for plots during each pulse. A minimum of 1 update is performed in each phase (2 for each pulse consisting of an illumination phase and a diffusion phase)
model.HS.mediaPropRecalcPeriod = 2; % Every N updates, the media properties will be recalculated (including, if needed, re-running MC and FMC steps)
model.HS.tempSensorPositions = [0 0.009 0.032
                                0 0.002 0.044
                                0 -0.007 0.059
                                0 -0.017 0.078]; % Each row is a temperature sensor's absolute [x y z] coordinates. Leave the matrix empty ([]) to disable temperature sensors.

% Execution, do not modify the next line:
model = simulateHeatDistribution(model);

plot(model,'HS');
plot(model,'MC');

%% Geometry function(s)
% A geometry function takes as input X,Y,Z matrices as returned by the
% "ndgrid" MATLAB function as well as any parameters the user may have
% provided in the definition of Ginput. It returns the media matrix M,
% containing numerical values indicating the media type (as defined in
% mediaPropertiesFunc) at each voxel location.
function M = geometryDefinition(X,Y,Z,parameters)
    M = ones(size(X)); % fill background with water
    M(0.2*Y + (Z - 0.03) > 0) = 2; % Thermochromic material
    M(-0.2*Y + (Z - 0.05) > 0) = 1; % Water
    M(Z > 0.07) = 2;
    M(Z > 0.09) = 1;
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
    mediaProperties(j).name  = 'water';
    mediaProperties(j).mua   = 0.00036; % [cm^-1]
    mediaProperties(j).mus   = 10; % [cm^-1]
    mediaProperties(j).g     = 1.0;
    mediaProperties(j).VHC   = 4.19; % [J cm^-3 K^-1]
    mediaProperties(j).TC    = 5.8e-3; % [W cm^-1 K^-1]
    
    j=2;
    mediaProperties(j).name  = 'thermochromic material';
    mediaProperties(j).mua = '150*max(1e-8,(40-T)/(40-20))'; % [cm^-1]
    mediaProperties(j).mus = 30; % [cm^-1]
    mediaProperties(j).g   = 0.8;
    mediaProperties(j).VHC = 2000; % [J cm^-3 K^-1]
    mediaProperties(j).TC  = 5e-3; % [W cm^-1 K^-1]
    mediaProperties(j).nBins = 150;
end