%% Decription
% This is a test to see if enabling all possible FR, FD and T dependences
% works.

%% Geometry definition
model = MCmatlab.model;

model.G.nx                = 100; % Number of bins in the x direction
model.G.ny                = 100; % Number of bins in the y direction
model.G.nz                = 100; % Number of bins in the z direction
model.G.Lx                = .1; % [cm] x size of simulation cuboid
model.G.Ly                = .1; % [cm] y size of simulation cuboid
model.G.Lz                = .1; % [cm] z size of simulation cuboid

model.G.mediaPropertiesFunc = @mediaPropertiesFunc; % Media properties defined as a function at the end of this file
model.G.geomFunc          = @geometryDefinition_AllFRTDependences; % Function to use for defining the distribution of media in the cuboid. Defined at the end of this m file.

plotMCmatlabGeom(model);

%% Monte Carlo simulation
model = reset(model,'MC'); % Only necessary if you want to run this section repeatedly, re-using previous G data

model.MC.simulationTimeRequested  = .05; % [min] Time duration of the simulation

model.MC.matchedInterfaces        = true; % Assumes all refractive indices are the same
model.MC.boundaryType             = 1; % 0: No escaping boundaries, 1: All cuboid boundaries are escaping, 2: Top cuboid boundary only is escaping
model.MC.wavelength               = 532; % [nm] Excitation wavelength, used for determination of optical properties for excitation light

model.MC.beam.beamType            = 5; % 0: Pencil beam, 1: Isotropically emitting point source, 2: Infinite plane wave, 3: Laguerre-Gaussian LG01 beam, 4: Radial-factorizable beam (e.g., a Gaussian beam), 5: X/Y factorizable beam (e.g., a rectangular LED emitter)
model.MC.beam.NF.XDistr           = 0; % X near field distribution - 0: Top-hat, 1: Gaussian, Array: Custom. Doesn't need to be normalized.
model.MC.beam.NF.XWidth           = .02; % [cm] X near field 1/e^2 radius if top-hat or Gaussian or half-width of the full distribution if custom
model.MC.beam.NF.YDistr           = 0; % Y near field distribution - 0: Top-hat, 1: Gaussian, Array: Custom. Doesn't need to be normalized.
model.MC.beam.NF.YWidth           = .01; % [cm] Y near field 1/e^2 radius if top-hat or Gaussian or half-width of the full distribution if custom
model.MC.beam.FF.XDistr           = 2; % X far field distribution - 0: Top-hat, 1: Gaussian, 2: Cosine (Lambertian), Array: Custom. Doesn't need to be normalized.
model.MC.beam.FF.XWidth           = pi/8; % [rad] X far field 1/e^2 half-angle if top-hat or Gaussian or half-angle of the full distribution if custom
model.MC.beam.FF.YDistr           = 2; % Y far field distribution - 0: Top-hat, 1: Gaussian, 2: Cosine (Lambertian), Array: Custom. Doesn't need to be normalized.
model.MC.beam.FF.YWidth           = pi/8; % [rad] Y far field 1/e^2 half-angle if top-hat or Gaussian or half-angle of the full distribution if custom
model.MC.beam.psi                 = pi/4; % [rad] (Default: 0) Axial rotation angle of beam, relevant only for XY distributed beams

model.MC.beam.xFocus              = 0; % [cm] x position of focus
model.MC.beam.yFocus              = 0; % [cm] y position of focus
model.MC.beam.zFocus              = 0; % [cm] z position of focus
model.MC.beam.theta               = 0; % [rad] Polar angle of beam center axis
model.MC.beam.phi                 = 0; % [rad] Azimuthal angle of beam center axis

model.MC.P                   = 2; % [W] Incident pulse peak power (in case of infinite plane waves, only the power incident upon the cuboid's top surface)

model.HS.Tinitial            = 5*rand(size(model.G.M_raw)); % [deg C] Initial temperature

% Execution, do not modify the next line:
model = runMonteCarlo(model);

plotMCmatlab(model);

%% Fluorescence Monte Carlo
model = reset(model,'FMC'); % Only necessary if you want to run this section repeatedly, re-using previous G and MC data

model.FMC.simulationTimeRequested = .1; % [min] Time duration of the simulation
model.FMC.wavelength              = 600; % [nm] Fluorescence wavelength, used for determination of optical properties for fluorescence light

model = runMonteCarlo(model,'fluorescence');

plotMCmatlab(model,'fluorescence');

%% Heat simulation
model.HS.largeTimeSteps      = true; % (Default: false) If true, calculations will be faster, but some voxel temperatures may be slightly less precise. Test for yourself whether this precision is acceptable for your application.

model.HS.heatBoundaryType    = 1; % 0: Insulating boundaries, 1: Constant-temperature boundaries (heat-sinked)
model.HS.durationOn          = 0.006; % [s] Pulse on-duration
model.HS.durationOff         = 0.004; % [s] Pulse off-duration
model.HS.durationEnd         = 0.003; % [s] Non-illuminated relaxation time to add to the end of the simulation to let temperature diffuse after the pulse train

model.HS.nPulses             = 2;
model.HS.plotTempLimits      = [0 45]; % [deg C] Expected range of temperatures, used only for setting the color scale in the plot
model.HS.nUpdates            = 10; % Number of times data is extracted for plots during each pulse. A minimum of 1 update is performed in each phase (2 for each pulse consisting of an illumination phase and a diffusion phase)
model.HS.mediaPropRecalcPeriod = 10; % Every N updates, the media properties will be recalculated (including, if needed, re-running MC and FMC steps)
model.HS.slicePositions      = [.5 0.6 1]; % Relative slice positions [x y z] for the 3D plots on a scale from 0 to 1
model.HS.tempSensorPositions = [0 0 0.005
                                0 0 0.015
                                0 0 0.035
                                0 0 0.065]; % Each row is a temperature sensor's absolute [x y z] coordinates. Leave the matrix empty ([]) to disable temperature sensors.

% Execution, do not modify the next line:
model = simulateHeatDistribution(model);

plotMCmatlabHeat(model);
plotMCmatlab(model);
plotMCmatlab(model,'fluorescence');

%% Post-processing

%% Geometry function(s)
% A geometry function takes as input X,Y,Z matrices as returned by the
% "ndgrid" MATLAB function as well as any parameters the user may have
% provided in the definition of Ginput. It returns the media matrix M,
% containing numerical values indicating the media type (as defined in
% mediaPropertiesFunc) at each voxel location.
function M = geometryDefinition_AllFRTDependences(X,Y,Z,parameters)
M = ones(size(X));
M(Z > 0.01) = 2;
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
mediaProperties(j).name  = 'material 1';
mediaProperties(j).mua   = '0.00036 + (T+FR)/100'; % [cm^-1]
mediaProperties(j).mus   = '10+FR/1000-T/100'; % [cm^-1]
mediaProperties(j).g     = '1./exp(T+FR/500)';
mediaProperties(j).PY    = '0.5/(1+FR/700+T/10)';
mediaProperties(j).VHC   = '4.19 + T/200';
mediaProperties(j).TC    = '5.8e-3 + T/25000';
mediaProperties(j).nBins = 20;

j=2;
mediaProperties(j).name  = 'material 2';
mediaProperties(j).mua   = '0.0073 + (T+2*FR)/60'; % [cm^-1]
mediaProperties(j).mus   = '50+FR/7000+T/150'; % [cm^-1]
mediaProperties(j).g     = '0.2*(0.5+0.5*tanh(T))+0.05*(0.5+0.5*tanh(FR/500))';
mediaProperties(j).QY    = '0.35-0.25*tanh((FR-250)/100)+0.1*tanh((T-10)/10)';
mediaProperties(j).VHC   = '2.19 + T/400';
mediaProperties(j).TC    = '10e-3 + T/150';
mediaProperties(j).nBins = 30;
end