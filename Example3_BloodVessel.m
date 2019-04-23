addpath([fileparts(matlab.desktop.editor.getActiveFilename) '/helperfuncs']); % The helperfuncs folder is added to the path for the duration of this MATLAB session

%% Decription
% This example simulates a collimated top hat beam of radius 300 µm
% incident on skin, with some gel (water) on the top. This example is
% constructed identically to that on the mcxyz website, except that photons
% escape on all boundaries and the voxel grid is only 100x100x100:
% https://omlc.org/software/mc/mcxyz/
%
% The found absorption distribution is then passed into the heat simulator,
% assuming the light is on for 5 pulses of 1 ms on time and 4 ms off time
% each, with 4 W of peak power. Some demonstration values of the Arrhenius
% E and A parameters for blood coagulation exist in getMediaProperties and
% are used to calculate the distribution of coagulated blood. Temperature
% sensors outputs and movie generation is also demonstrated.

%% Geometry definition
clear Ginput
Ginput.matchedInterfaces = true; % Assumes all refractive indices are 1
Ginput.boundaryType      = 1; % 0: No escaping boundaries, 1: All cuboid boundaries are escaping, 2: Top cuboid boundary only is escaping

Ginput.wavelength        = 532; % [nm] Excitation wavelength, used for determination of optical properties for excitation light

Ginput.nx                = 200; % Number of bins in the x direction
Ginput.ny                = 200; % Number of bins in the y direction
Ginput.nz                = 200; % Number of bins in the z direction
Ginput.Lx                = .1; % [cm] x size of simulation cuboid
Ginput.Ly                = .1; % [cm] y size of simulation cuboid
Ginput.Lz                = .1; % [cm] z size of simulation cuboid

Ginput.GeomFunc          = @GeometryDefinition_BloodVessel; % Function to use for defining the distribution of media in the cuboid. Defined at the end of this m file.

% Execution, do not modify the next two lines:
Goutput = defineGeometry(Ginput);
plotMCmatlabGeom(Goutput);

%% Monte Carlo simulation
clear MCinput
MCinput.simulationTime           = .1; % [min] Time duration of the simulation

MCinput.Beam.beamType            = 2; % 0: Pencil beam, 1: Isotropically emitting point source, 2: Infinite plane wave, 3: Gaussian focus, Gaussian far field beam, 4: Gaussian focus, top-hat far field beam, 5: Top-hat focus, Gaussian far field beam, 6: Top-hat focus, top-hat far field beam, 7: Laguerre-Gaussian LG01 beam
MCinput.Beam.xFocus              = 0; % [cm] x position of focus
MCinput.Beam.yFocus              = 0; % [cm] y position of focus
MCinput.Beam.zFocus              = 0; % [cm] z position of focus
MCinput.Beam.theta               = 0; % [rad] Polar angle of beam center axis
MCinput.Beam.phi                 = 0; % [rad] Azimuthal angle of beam center axis
MCinput.Beam.waist               = 0.03; % [cm] Beam waist 1/e^2 radius
MCinput.Beam.divergence          = 0/180*pi; % [rad] Beam divergence 1/e^2 half-angle of beam (for a diffraction limited Gaussian beam, this is G.wavelength*1e-9/(pi*MCinput.Beam.waist*1e-2))

% Execution, do not modify the next three lines:
MCinput.G = Goutput;
MCoutput = runMonteCarlo(MCinput);
plotMCmatlab(MCinput,MCoutput);

%% Heat simulation
HSinput.useAllCPUs          = true; % If false, MCmatlab will leave one processor unused. Useful for doing other work on the PC while simulations are running.
HSinput.makeMovie           = false; % Requires silentMode = false.

HSinput.heatBoundaryType    = 0; % 0: Insulating boundaries, 1: Constant-temperature boundaries (heat-sinked)
HSinput.P                   = 4; % [W] Incident pulse peak power (in case of infinite plane waves, only the power incident upon the cuboid's top surface)
HSinput.durationOn          = 0.0; % [s] Pulse on-duration
HSinput.durationOff         = 0.01; % [s] Pulse off-duration
HSinput.durationEnd         = 0.0; % [s] Non-illuminated relaxation time to add to the end of the simulation to let temperature diffuse after the pulse train
HSinput.initialTemp         = 50; % [deg C] Initial temperature

HSinput.nPulses             = 1; % Number of consecutive pulses, each with an illumination phase and a diffusion phase. If simulating only illumination or only diffusion, use n_pulses = 1.

HSinput.plotTempLimits      = [1 200]; % [deg C] Expected range of temperatures, used only for setting the color scale in the plot
HSinput.nUpdates            = 1; % Number of times data is extracted for plots during each pulse. A minimum of 1 update is performed in each phase (2 for each pulse consisting of an illumination phase and a diffusion phase)
HSinput.slicePositions      = [.5 0.6 1]; % Relative slice positions [x y z] for the 3D plots on a scale from 0 to 1
HSinput.tempSensorPositions = [0 0 0.038
                               0 0 0.04
                               0 0 0.042
                               0 0 0.044]; % Each row is a temperature sensor's absolute [x y z] coordinates. Leave the matrix empty ([]) to disable temperature sensors.

% Execution, do not modify the next four lines:
HSinput.G = Goutput;
HSinput.MCoutput = MCoutput;
HSoutput = simulateHeatDistribution(HSinput);
plotMCmatlabHeat(HSinput,HSoutput);

max(HSoutput.Omega(:))

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
M = 5*ones(size(X)); % fill background with water (gel)
% M = 2*ones(size(X)); % fill background with water (gel)
% M(Z > zsurf) = 4; % epidermis
% M(Z > zsurf + epd_thick) = 5; % dermis
% M(X.^2 + (Z - (zsurf + vesseldepth)).^2 < vesselradius^2) = 6; % blood
end
