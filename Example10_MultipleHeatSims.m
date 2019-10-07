addpath([fileparts(matlab.desktop.editor.getActiveFilename) '/helperfuncs']); % The helperfuncs folder is added to the path for the duration of this MATLAB session

%% Decription
% This example is very similar to example 3, the geometry and light
% propagation part being the same. It illustrates the ability to sequence
% together multiple pulse trains in the form of separate heat simulations,
% carrying over the temperature and thermal damage from one simulation to
% the next. A movie is generated that shows all the heat simulations.
% 
% The first heat simulation has a single light pulse of 2 ms on-time
% followed by 3 ms off-time. The second heat simulation adds nine more
% pulses onto that with 0.5 ms on-time and 4.5 ms off-time.
% 
% In these heat simulations, largeTimeSteps is set to true and nUpdates set
% to 10 rather than 100, which speeds up the simulations considerably.

%% Geometry definition
clear Ginput
Ginput.nx                = 100; % Number of bins in the x direction
Ginput.ny                = 100; % Number of bins in the y direction
Ginput.nz                = 100; % Number of bins in the z direction
Ginput.Lx                = .1; % [cm] x size of simulation cuboid
Ginput.Ly                = .1; % [cm] y size of simulation cuboid
Ginput.Lz                = .1; % [cm] z size of simulation cuboid

Ginput.mediaPropertiesFunc = @mediaPropertiesFunc; % Media properties defined as a function at the end of this file
Ginput.GeomFunc          = @GeometryDefinition_BloodVessel; % Function to use for defining the distribution of media in the cuboid. Defined at the end of this m file.

% Execution, do not modify the next line:
Goutput = defineGeometry(Ginput);

plotMCmatlabGeom(Goutput);

%% Monte Carlo simulation
clear MCinput
MCinput.simulationTime           = .1; % [min] Time duration of the simulation

MCinput.matchedInterfaces        = true; % Assumes all refractive indices are 1
MCinput.boundaryType             = 1; % 0: No escaping boundaries, 1: All cuboid boundaries are escaping, 2: Top cuboid boundary only is escaping
MCinput.wavelength               = 532; % [nm] Excitation wavelength, used for determination of optical properties for excitation light

MCinput.Beam.beamType            = 6; % 0: Pencil beam, 1: Isotropically emitting point source, 2: Infinite plane wave, 3: Gaussian focus, Gaussian far field beam, 4: Gaussian focus, top-hat far field beam, 5: Top-hat focus, Gaussian far field beam, 6: Top-hat focus, top-hat far field beam, 7: Laguerre-Gaussian LG01 beam
MCinput.Beam.xFocus              = 0; % [cm] x position of focus
MCinput.Beam.yFocus              = 0; % [cm] y position of focus
MCinput.Beam.zFocus              = 0; % [cm] z position of focus
MCinput.Beam.theta               = 0; % [rad] Polar angle of beam center axis
MCinput.Beam.phi                 = 0; % [rad] Azimuthal angle of beam center axis
MCinput.Beam.waist               = 0.03; % [cm] Beam waist 1/e^2 radius
MCinput.Beam.divergence          = 0/180*pi; % [rad] Beam divergence 1/e^2 half-angle of beam (for a diffraction limited Gaussian beam, this is G.wavelength*1e-9/(pi*MCinput.Beam.waist*1e-2))

% Execution, do not modify the next two lines:
MCinput.G = Goutput;
MCoutput = runMonteCarlo(MCinput);

% plotMCmatlab(MCinput,MCoutput);

%% First pulse heat simulation
clear HSinput
HSinput.useAllCPUs          = true; % If false, MCmatlab will leave one processor unused. Useful for doing other work on the PC while simulations are running.
HSinput.makeMovie           = true; % Requires silentMode = false.
HSinput.deferMovieWrite     = true; % (Default: false) If true, will not write the movie to file upon completion, but will store the raw frames in HSoutput for future writing
HSinput.largeTimeSteps      = true; % (Default: false) If true, calculations will be faster, but some voxel temperatures may be slightly less precise. Test for yourself whether this precision is acceptable for your application.

HSinput.heatBoundaryType    = 0; % 0: Insulating boundaries, 1: Constant-temperature boundaries (heat-sinked)
HSinput.P                   = 4; % [W] Incident pulse peak power (in case of infinite plane waves, only the power incident upon the cuboid's top surface)
HSinput.durationOn          = 0.002; % [s] Pulse on-duration
HSinput.durationOff         = 0.003; % [s] Pulse off-duration
HSinput.durationEnd         = 0.000; % [s] Non-illuminated relaxation time to add to the end of the simulation to let temperature diffuse after the pulse train
HSinput.initialTemp         = 37; % [deg C] Initial temperature

HSinput.nPulses             = 1; % Number of consecutive pulses, each with an illumination phase and a diffusion phase. If simulating only illumination or only diffusion, use n_pulses = 1.

HSinput.plotTempLimits      = [37 100]; % [deg C] Expected range of temperatures, used only for setting the color scale in the plot
HSinput.nUpdates            = 10; % Number of times data is extracted for plots during each pulse. A minimum of 1 update is performed in each phase (2 for each pulse consisting of an illumination phase and a diffusion phase)
HSinput.slicePositions      = [.5 0.6 1]; % Relative slice positions [x y z] for the 3D plots on a scale from 0 to 1
HSinput.tempSensorPositions = [0 0 0.038
                               0 0 0.04
                               0 0 0.042
                               0 0 0.044]; % Each row is a temperature sensor's absolute [x y z] coordinates. Leave the matrix empty ([]) to disable temperature sensors.

% Execution, do not modify the next three lines:
HSinput.G = Goutput;
HSinput.MCoutput = MCoutput;
HSoutput = simulateHeatDistribution(HSinput);

% plotMCmatlabHeat(HSinput,HSoutput);

%% Remaining pulse train heat simulation
% The HSinput struct carries over from the previous section and is modified slightly
HSinput.prevHSoutput        = HSoutput; % The HSoutput of the previous heat simulation, to use as the new initial condition
HSinput.deferMovieWrite     = false; % (Default: false) If true, will not write the movie to file upon completion, but will store the raw frames in HSoutput for future writing
HSinput.durationOn          = 0.0005; % [s] Pulse on-duration
HSinput.durationOff         = 0.0045; % [s] Pulse off-duration
HSinput.durationEnd         = 0.01; % [s] Non-illuminated relaxation time to add to the end of the simulation to let temperature diffuse after the pulse train

HSinput.nPulses             = 9; % Number of consecutive pulses, each with an illumination phase and a diffusion phase. If simulating only illumination or only diffusion, use n_pulses = 1.

% Execution, do not modify the next three lines:
HSinput.G = Goutput;
HSinput.MCoutput = MCoutput;
HSoutput = simulateHeatDistribution(HSinput);

plotMCmatlabHeat(HSinput,HSoutput);

%% Post-processing

%% Geometry function(s)
% A geometry function takes as input X,Y,Z matrices as returned by the
% "ndgrid" MATLAB function as well as any parameters the user may have
% provided in the definition of Ginput. It returns the media matrix M,
% containing numerical values indicating the media type (as defined in
% getMediaProperties) at each voxel location.
function M = GeometryDefinition_BloodVessel(X,Y,Z,parameters)
% Blood vessel example:
zsurf = 0.01;
epd_thick = 0.006;
vesselradius  = 0.0100;
vesseldepth = 0.04;
M = ones(size(X)); % fill background with water (gel)
M(Z > zsurf) = 2; % epidermis
M(Z > zsurf + epd_thick) = 3; % dermis
M(X.^2 + (Z - (zsurf + vesseldepth)).^2 < vesselradius^2) = 4; % blood
end

%% Media Properties function
% The media properties function defines all the optical and thermal
% properties of the media involved by constructing and returning a
% "mediaProperties" struct with various fields. As its input, the function
% takes the wavelength as well as any other parameters you might specify
% above in the model file, for example parameters that you might loop over
% in a for loop.
function mediaProperties = mediaPropertiesFunc(wavelength,parameters)
load spectralLIB.mat
MU(:,1) = interp1(nmLIB,muaoxy,wavelength);
MU(:,2) = interp1(nmLIB,muadeoxy,wavelength);
MU(:,3) = interp1(nmLIB,muawater,wavelength);
MU(:,4) = interp1(nmLIB,muamel,wavelength);

j=1;
mediaProperties(j).name  = 'water';
mediaProperties(j).mua   = 0.00036;
mediaProperties(j).mus   = 10;
mediaProperties(j).g     = 1.0;
mediaProperties(j).n     = 1.3;
mediaProperties(j).VHC   = 4.19;
mediaProperties(j).TC    = 5.8e-3;

j=2;
mediaProperties(j).name  = 'epidermis';
B = 0;
S = 0.75;
W = 0.75;
Me = 0.03;
musp500 = 40;
fray    = 0.0;
bmie    = 1.0;
gg      = 0.90;
musp = musp500*(fray*(wavelength/500).^-4 + (1-fray)*(wavelength/500).^-bmie);
X = [B*S B*(1-S) W Me]';
mediaProperties(j).mua = MU*X;
mediaProperties(j).mus = musp/(1-gg);
mediaProperties(j).g   = gg;
mediaProperties(j).n   = 1.3;
mediaProperties(j).VHC = 3391*1.109e-3;
mediaProperties(j).TC  = 0.37e-2;

j=3;
mediaProperties(j).name = 'dermis';
B = 0.002;
S = 0.67;
W = 0.65;
Me = 0;
musp500 = 42.4;
fray    = 0.62;
bmie    = 1.0;
gg      = 0.90;
musp = musp500*(fray*(wavelength/500).^-4 + (1-fray)*(wavelength/500).^-bmie);
X = [B*S B*(1-S) W Me]';
mediaProperties(j).mua = MU*X;
mediaProperties(j).mus = musp/(1-gg);
mediaProperties(j).g   = gg;
mediaProperties(j).n   = 1.3;
mediaProperties(j).VHC = 3391*1.109e-3;
mediaProperties(j).TC  = 0.37e-2;

j=4;
mediaProperties(j).name  = 'blood';
B       = 1.00;
S       = 0.75;
W       = 0.95;
Me      = 0;
musp500 = 10;
fray    = 0.0;
bmie    = 1.0;
gg      = 0.90;
musp = musp500*(fray*(wavelength/500).^-4 + (1-fray)*(wavelength/500).^-bmie);
X = [B*S B*(1-S) W Me]';
mediaProperties(j).mua = MU*X;
mediaProperties(j).mus = musp/(1-gg);
mediaProperties(j).g   = gg;
mediaProperties(j).n   = 1.3;
mediaProperties(j).VHC = 3617*1.050e-3;
mediaProperties(j).TC  = 0.52e-2;
mediaProperties(j).E   = 422.5e3; % J/mol    PLACEHOLDER DATA ONLY
mediaProperties(j).A   = 7.6e66; % 1/s        PLACEHOLDER DATA ONLY
end