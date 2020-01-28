function model = clearMCmatlabModel(model,string)
% NaN means not yet specified or calculated and no default value exists

if strcmp(string,'G')
  %% Geometry parameters that can be specified by the user
  model.G.silentMode = false;
  
  model.G.nx = NaN;
  model.G.ny = NaN;
  model.G.nz = NaN;
  model.G.Lx = NaN;
  model.G.Ly = NaN;
  model.G.Lz = NaN;
  
  model.G.mediaPropertiesFunc = NaN;
  model.G.geomFunc = NaN;
  model.G.mediaPropParams = {};
  model.G.geomFuncParams = {};
  
  %% Geometry parameters that are calculated
  model.G.M_raw = NaN;
  
  model.G.dx = NaN;
  model.G.dy = NaN;
  model.G.dz = NaN;
  model.G.x = NaN;
  model.G.y = NaN;
  model.G.z = NaN;
end

if strcmp(string,'MC')
  %% Monte Carlo parameters that can be specified by the user
  model.MC.useGPU = false;
%   model.MC.useGPU = true;

  model.MC.simulationTimeRequested = 0.1;
  model.MC.nPhotonsRequested = NaN;
  model.MC.silentMode = false;
  model.MC.useAllCPUs = false;
  model.MC.calcNFR = true;
  model.MC.calcNFRdet = false;
  model.MC.nExamplePaths = 0;
  model.MC.farFieldRes = 0;
  
  model.MC.matchedInterfaces = true;
  model.MC.boundaryType = 1;
  model.MC.wavelength = NaN;
  
  model.MC.P = NaN;
  model.MC.FRinitial = NaN;
  model.MC.FR = NaN;
  model.MC.FRdepIterations = 20;
  
  % Beam parameters
  model.MC.beam.beamType = NaN;
  model.MC.beam.xFocus = 0;
  model.MC.beam.yFocus = 0;
  model.MC.beam.zFocus= 0;
  model.MC.beam.theta = 0;
  model.MC.beam.phi = 0;
  model.MC.beam.psi = 0;

  model.MC.beam.NF.radialDistr = NaN;
  model.MC.beam.NF.radialWidth = NaN;
  model.MC.beam.NF.XDistr = NaN;
  model.MC.beam.NF.XWidth = NaN;
  model.MC.beam.NF.YDistr = NaN;
  model.MC.beam.NF.YWidth = NaN;

  model.MC.beam.FF.radialDistr = NaN;
  model.MC.beam.FF.radialWidth = NaN;
  model.MC.beam.FF.XDistr = NaN;
  model.MC.beam.FF.XWidth = NaN;
  model.MC.beam.FF.YDistr = NaN;
  model.MC.beam.FF.YWidth = NaN;

  % Light collector parameters
  model.MC.useLightCollector = false;
  
  model.MC.LC.x = NaN;
  model.MC.LC.y = NaN;
  model.MC.LC.z = NaN;
  model.MC.LC.theta = NaN;
  model.MC.LC.phi = NaN;
  
  model.MC.LC.f = NaN;
  model.MC.LC.diam = NaN;
  model.MC.LC.fieldSize = NaN;
  model.MC.LC.NA = NaN;
  model.MC.LC.res = NaN;
  
  model.MC.LC.tStart = NaN;
  model.MC.LC.tEnd = NaN;
  model.MC.LC.nTimeBins = 0;
  
  %% Monte Carlo parameters that are calculated
  model.MC.simulationTime = NaN;
  model.MC.nPhotons = NaN;
  model.MC.nThreads = NaN;
  
  model.MC.mediaProperties_funcHandles = NaN; % Wavelength-dependent
  model.MC.mediaProperties = NaN; % Wavelength- and splitting-dependent
  model.MC.FRdependent = NaN;
  model.MC.FDdependent = NaN;
  model.MC.Tdependent = NaN;
  model.MC.M = NaN; % Splitting-dependent
  model.MC.RI = NaN;
  
  model.MC.examplePaths = NaN;
  
  model.MC.NFR = NaN; % Normalized Fluence Rate
  model.MC.NFRdet = NaN;
  
  model.MC.farField = NaN;
  model.MC.farFieldTheta = NaN;
  model.MC.farFieldPhi = NaN;
  
  model.MC.LC.image = NaN;
  model.MC.LC.X = NaN;
  model.MC.LC.Y = NaN;
  
  model.MC.NI_xpos = NaN; % Normalized irradiance on the boundary in the positive x direction
  model.MC.NI_xneg = NaN;
  model.MC.NI_ypos = NaN;
  model.MC.NI_yneg = NaN;
  model.MC.NI_zpos = NaN;
  model.MC.NI_zneg = NaN;
end

if strcmp(string,'FMC')
  %% Fluorescence Monte Carlo parameters that can be specified by the user
  model.FMC.useGPU = false;
%   model.FMC.useGPU = true;

  model.FMC.simulationTimeRequested = 0.1;
  model.FMC.nPhotonsRequested = NaN;
  model.FMC.silentMode = false;
  model.FMC.useAllCPUs = false;
  model.FMC.calcNFR = true;
  model.FMC.calcNFRdet = false;
  model.FMC.nExamplePaths = 0;
  model.FMC.farFieldRes = 0;
  
  model.FMC.matchedInterfaces = true;
  model.FMC.boundaryType = 1;
  model.FMC.wavelength = NaN;
  
  % Light collector parameters
  model.FMC.useLightCollector = false;
  
  model.FMC.LC.x = NaN;
  model.FMC.LC.y = NaN;
  model.FMC.LC.z = NaN;
  model.FMC.LC.theta = NaN;
  model.FMC.LC.phi = NaN;
  
  model.FMC.LC.f = NaN;
  model.FMC.LC.diam = NaN;
  model.FMC.LC.fieldSize = NaN;
  model.FMC.LC.NA = NaN;
  model.FMC.LC.res = NaN;
  
  %% Fluorescence Monte Carlo parameters that are calculated
  model.FMC.simulationTime = NaN;
  model.FMC.nPhotons = NaN;
  model.FMC.nThreads = NaN;
  
  model.FMC.mediaProperties_funcHandles = NaN; % Wavelength-dependent
  model.FMC.mediaProperties = NaN; % Wavelength- and splitting-dependent
  model.FMC.FRdependent = NaN;
  model.FMC.FDdependent = NaN;
  model.FMC.Tdependent = NaN;
  model.FMC.M = NaN; % Splitting-dependent
  model.FMC.RI = NaN;
  
  model.FMC.examplePaths = NaN;
  
  model.FMC.NFR = NaN; % Normalized Fluence Rate
  model.FMC.NFRdet = NaN;
  
  model.FMC.farField = NaN;
  model.FMC.farFieldTheta = NaN;
  model.FMC.farFieldPhi = NaN;
  
  model.FMC.beam.sourceDistribution = NaN;
  
  model.FMC.LC.image = NaN;
  model.FMC.LC.X = NaN;
  model.FMC.LC.Y = NaN;
  
  model.FMC.NI_xpos = NaN; % Normalized irradiance on the boundary in the positive x direction
  model.FMC.NI_xneg = NaN;
  model.FMC.NI_ypos = NaN;
  model.FMC.NI_yneg = NaN;
  model.FMC.NI_zpos = NaN;
  model.FMC.NI_zneg = NaN;
end

if strcmp(string,'HS')
  %% Heat simulation parameters that can be specified by the user
  model.HS.silentMode = false;
  model.HS.useAllCPUs = false;
  model.HS.useGPU     = false;
  model.HS.makeMovie = false;
  model.HS.largeTimeSteps = false;
  model.HS.deferMovieWrite = false;
  
  model.HS.mediaProperties_funcHandles = NaN;
  model.HS.mediaProperties = NaN; % Splitting-dependent
  
  model.HS.heatBoundaryType = 0; % 0: Insulating boundaries, 1: Constant-temperature boundaries (heat-sinked)
  model.HS.durationOn = NaN; % [s] Pulse on-duration
  model.HS.durationOff = NaN; % [s] Pulse off-duration
  model.HS.durationEnd = NaN; % [s] Non-illuminated relaxation time to add to the end of the simulation to let temperature diffuse after the pulse train
  
  model.HS.nPulses = 1; % Number of consecutive pulses, each with an illumination phase and a diffusion phase. If simulating only illumination or only diffusion, use n_pulses = 1.
  
  model.HS.plotTempLimits = NaN; % [deg C] Expected range of temperatures, used only for setting the color scale in the plot
  model.HS.nUpdates = 10; % Number of times data is extracted for plots during each pulse. A minimum of 1 update is performed in each phase (2 for each pulse consisting of an illumination phase and a diffusion phase)
  model.HS.mediaPropRecalcPeriod = 1; % If necessary, recalculates the fluence rate and thus any fluence rate dependent optical and thermal properties every N updates
  
  model.HS.slicePositions = [0.5 1 1];
  model.HS.tempSensorPositions = [];
  model.HS.Tinitial = NaN;
  
  %% Heat simulation parameters that are calculated
  model.HS.Tdependent = NaN;
  model.HS.FDdependent = NaN;

  model.HS.T = NaN; % Final temperature
  model.HS.Omega = single(NaN);
  model.HS.maxMediaTemps = NaN;
  model.HS.sensorsTimeVector = NaN;
  model.HS.sensorTemps = NaN;
  model.HS.movieFrames = NaN;
end

end