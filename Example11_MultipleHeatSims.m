%% Decription
% This example is very similar to example 4, the geometry and light
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
MCmatlab.closeMCmatlabFigures();
model = MCmatlab.model;

model.G.nx                = 100; % Number of bins in the x direction
model.G.ny                = 100; % Number of bins in the y direction
model.G.nz                = 100; % Number of bins in the z direction
model.G.Lx                = .1; % [cm] x size of simulation cuboid
model.G.Ly                = .1; % [cm] y size of simulation cuboid
model.G.Lz                = .1; % [cm] z size of simulation cuboid

model.G.mediaPropertiesFunc = @mediaPropertiesFunc; % Media properties defined as a function at the end of this file
model.G.geomFunc          = @geometryDefinition; % Function to use for defining the distribution of media in the cuboid. Defined at the end of this m file.

model = plot(model,'G');

%% Monte Carlo simulation
model.MC.simulationTimeRequested  = .1; % [min] Time duration of the simulation

model.MC.matchedInterfaces        = true; % Assumes all refractive indices are the same
model.MC.boundaryType             = 1; % 0: No escaping boundaries, 1: All cuboid boundaries are escaping, 2: Top cuboid boundary only is escaping, 3: Top and bottom boundaries are escaping, while the side boundaries are cyclic
model.MC.wavelength               = 532; % [nm] Excitation wavelength, used for determination of optical properties for excitation light

model.MC.lightSource.sourceType   = 4; % 0: Pencil beam, 1: Isotropically emitting line or point source, 2: Infinite plane wave, 3: Laguerre-Gaussian LG01 beam, 4: Radial-factorizable beam (e.g., a Gaussian beam), 5: X/Y factorizable beam (e.g., a rectangular LED emitter)
model.MC.lightSource.focalPlaneIntensityDistribution.radialDistr = 0; % Radial focal plane intensity distribution - 0: Top-hat, 1: Gaussian, Array: Custom. Doesn't need to be normalized.
model.MC.lightSource.focalPlaneIntensityDistribution.radialWidth = .03; % [cm] Radial focal plane 1/e^2 radius if top-hat or Gaussian or half-width of the full distribution if custom
model.MC.lightSource.angularIntensityDistribution.radialDistr = 0; % Radial angular intensity distribution - 0: Top-hat, 1: Gaussian, 2: Cosine (Lambertian), Array: Custom. Doesn't need to be normalized.
model.MC.lightSource.angularIntensityDistribution.radialWidth = 0; % [rad] Radial angular 1/e^2 half-angle if top-hat or Gaussian or half-angle of the full distribution if custom. For a diffraction limited Gaussian beam, this should be set to model.MC.wavelength*1e-9/(pi*model.MC.lightSource.focalPlaneIntensityDistribution.radialWidth*1e-2))
model.MC.lightSource.xFocus       = 0; % [cm] x position of focus
model.MC.lightSource.yFocus       = 0; % [cm] y position of focus
model.MC.lightSource.zFocus       = 0; % [cm] z position of focus
model.MC.lightSource.theta        = 0; % [rad] Polar angle of beam center axis
model.MC.lightSource.phi          = 0; % [rad] Azimuthal angle of beam center axis


model = runMonteCarlo(model);
model = plot(model,'MC');

%% First pulse heat simulation
model.MC.P                   = 4; % [W] Incident pulse peak power (in case of infinite plane waves, only the power incident upon the cuboid's top surface)

model.HS.useAllCPUs          = true; % If false, MCmatlab will leave one processor unused. Useful for doing other work on the PC while simulations are running.
model.HS.makeMovie           = true; % Requires silentMode = false.
model.HS.deferMovieWrite     = true; % (Default: false) If true, will not write the movie to file upon completion, but will store the raw frames in HSoutput for future writing
model.HS.largeTimeSteps      = true; % (Default: false) If true, calculations will be faster, but some voxel temperatures may be slightly less precise. Test for yourself whether this precision is acceptable for your application.

model.HS.heatBoundaryType    = 0; % 0: Insulating boundaries, 1: Constant-temperature boundaries (heat-sinked)
model.HS.durationOn          = 0.002; % [s] Pulse on-duration
model.HS.durationOff         = 0.003; % [s] Pulse off-duration
model.HS.durationEnd         = 0.000; % [s] Non-illuminated relaxation time to add to the end of the simulation to let temperature diffuse after the pulse train
model.HS.T                   = 37; % [deg C] Initial temperature

model.HS.nPulses             = 1; % Number of consecutive pulses, each with an illumination phase and a diffusion phase. If simulating only illumination or only diffusion, use nPulses = 1.

model.HS.plotTempLimits      = [37 100]; % [deg C] Expected range of temperatures, used only for setting the color scale in the plot
model.HS.nUpdates            = 10; % Number of times data is extracted for plots during each pulse. A minimum of 1 update is performed in each phase (2 for each pulse consisting of an illumination phase and a diffusion phase)
model.HS.slicePositions      = [.5 0.6 1]; % Relative slice positions [x y z] for the 3D plots on a scale from 0 to 1
model.HS.tempSensorPositions = [0 0 0.038
                                0 0 0.04
                                0 0 0.042
                                0 0 0.044]; % Each row is a temperature sensor's absolute [x y z] coordinates. Leave the matrix empty ([]) to disable temperature sensors.


model = simulateHeatDistribution(model);

%% Remaining pulse train heat simulation
% The model.HS struct is modified slightly but keeps its data from the
% previous pulse (HS.T, HS.sensorTemps, HS.Omega, etc.)
model.HS.deferMovieWrite     = false; % (Default: false) If true, will not write the movie to file upon completion, but will store the raw frames in HSoutput for future writing
model.HS.durationOn          = 0.0005; % [s] Pulse on-duration
model.HS.durationOff         = 0.0045; % [s] Pulse off-duration
model.HS.durationEnd         = 0.01; % [s] Non-illuminated relaxation time to add to the end of the simulation to let temperature diffuse after the pulse train

model.HS.nPulses             = 9; % Number of consecutive pulses, each with an illumination phase and a diffusion phase. If simulating only illumination or only diffusion, use nPulses = 1.


model = simulateHeatDistribution(model);
model = plot(model,'HS');

%% Geometry function(s) (see readme for details)
function M = geometryDefinition(X,Y,Z,parameters)
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

%% Media Properties function (see readme for details)
function mediaProperties = mediaPropertiesFunc(parameters)
  mediaProperties = MCmatlab.mediumProperties;

  j=1;
  mediaProperties(j).name  = 'water';
  mediaProperties(j).mua   = 0.00036; % [cm^-1]
  mediaProperties(j).mus   = 10; % [cm^-1]
  mediaProperties(j).g     = 1.0;
  mediaProperties(j).VHC   = 4.19; % [J cm^-3 K^-1]
  mediaProperties(j).TC    = 5.8e-3; % [W cm^-1 K^-1]

  j=2;
  mediaProperties(j).name  = 'epidermis';
  mediaProperties(j).mua = @func_mua2;
  function mua = func_mua2(wavelength)
    B = 0; % Blood content
    S = 0.75; % Blood oxygen saturation
    W = 0.75; % Water content
    M = 0.03; % Melanin content
    F = 0; % Fat content
    mua = calc_mua(wavelength,S,B,W,F,M); % Jacques "Optical properties of biological tissues: a review" eq. 12
  end

  mediaProperties(j).mus = @func_mus2;
  function mus = func_mus2(wavelength)
    aPrime = 40; % musPrime at 500 nm
    fRay = 0; % Fraction of scattering due to Rayleigh scattering
    bMie = 1; % Scattering power for Mie scattering
    g = 0.9; % Scattering anisotropy
    mus = calc_mus(wavelength,aPrime,fRay,bMie,g); % Jacques "Optical properties of biological tissues: a review" eq. 2
  end
  mediaProperties(j).g   = 0.9;
  mediaProperties(j).VHC = 3391*1.109e-3; % [J cm^-3 K^-1]
  mediaProperties(j).TC  = 0.37e-2; % [W cm^-1 K^-1]

  j=3;
  mediaProperties(j).name = 'dermis';
  mediaProperties(j).mua = @func_mua3;
  function mua = func_mua3(wavelength)
    B = 0.002; % Blood content
    S = 0.67; % Blood oxygen saturation
    W = 0.65; % Water content
    M = 0; % Melanin content
    F = 0; % Fat content
    mua = calc_mua(wavelength,S,B,W,F,M); % Jacques "Optical properties of biological tissues: a review" eq. 12
  end

  mediaProperties(j).mus = @func_mus3;
  function mus = func_mus3(wavelength)
    aPrime = 42.4; % musPrime at 500 nm
    fRay = 0.62; % Fraction of scattering due to Rayleigh scattering
    bMie = 1; % Scattering power for Mie scattering
    g = 0.9; % Scattering anisotropy
    mus = calc_mus(wavelength,aPrime,fRay,bMie,g); % Jacques "Optical properties of biological tissues: a review" eq. 2
  end
  mediaProperties(j).g   = 0.9;
  mediaProperties(j).VHC = 3391*1.109e-3; % [J cm^-3 K^-1]
  mediaProperties(j).TC  = 0.37e-2; % [W cm^-1 K^-1]

  j=4;
  mediaProperties(j).name  = 'blood';
  mediaProperties(j).mua = @func_mua4;
  function mua = func_mua4(wavelength)
    B = 1; % Blood content
    S = 0.75; % Blood oxygen saturation
    W = 0.95; % Water content
    M = 0; % Melanin content
    F = 0; % Fat content
    mua = calc_mua(wavelength,S,B,W,F,M); % Jacques "Optical properties of biological tissues: a review" eq. 12
  end
  
  mediaProperties(j).mus = @func_mus4;
  function mus = func_mus4(wavelength)
    aPrime = 10; % musPrime at 500 nm
    fRay = 0; % Fraction of scattering due to Rayleigh scattering
    bMie = 1; % Scattering power for Mie scattering
    g = 0.9; % Scattering anisotropy
    mus = calc_mus(wavelength,aPrime,fRay,bMie,g); % Jacques "Optical properties of biological tissues: a review" eq. 2
  end
  mediaProperties(j).g   = 0.9;
  mediaProperties(j).VHC = 3617*1.050e-3; % [J cm^-3 K^-1]
  mediaProperties(j).TC  = 0.52e-2; % [W cm^-1 K^-1]
  mediaProperties(j).E   = 422.5e3; % J/mol    PLACEHOLDER DATA ONLY
  mediaProperties(j).A   = 7.6e66; % 1/s        PLACEHOLDER DATA ONLY
end