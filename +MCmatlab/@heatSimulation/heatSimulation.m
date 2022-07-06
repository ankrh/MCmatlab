classdef heatSimulation
  %HEATSIMULATION  This class includes all properties and methods
  %related to the heat simulation in an MCmatlab.model.

  properties
    %% Input properties
    silentMode (1,1) logical = false % Disables command window text and progress indication
    useAllCPUs (1,1) logical = false % If false, MCmatlab will leave one processor unused. Useful for doing other work on the PC while simulations are running.
    useGPU (1,1) logical = false % Use CUDA acceleration for NVIDIA GPUs
    makeMovie (1,1) logical = false % Requires silentMode = false.
    deferMovieWrite (1,1) logical = false
    largeTimeSteps (1,1) logical = false % If true, calculations will be faster, but some voxel temperatures may be slightly less precise. Test for yourself whether this precision is acceptable for your application.

    heatBoundaryType (1,1) double {mustBeInteger, mustBeInRange(heatBoundaryType,0,1)} = 0 % 0: Insulating boundaries, 1: Constant-temperature boundaries (heat-sinked)
    durationOn (1,1) double {mustBeNonnegative} = 0 % [s] Pulse on-duration
    durationOff (1,1) double {mustBeNonnegative} = 0 % [s] Pulse off-duration
    durationEnd (1,1) double {mustBeNonnegative} = 0 % [s] Non-illuminated relaxation time to add to the end of the simulation to let temperature diffuse after the pulse train
    Tinitial (:,:,:) double {mustBeFiniteOrNaNScalar} = NaN % [deg C] Initial temperature, can be scalar or 3D array

    nPulses (1,1) double {mustBeInteger, mustBePositive} = 1 % Number of consecutive pulses, each with an illumination phase and a diffusion phase. If simulating only illumination or only diffusion, use n_pulses = 1.

    plotTempLimits (1,2) {mustBeAllFiniteOrAllNaN} = [NaN NaN] % [deg C] Expected range of temperatures, used only for setting the color scale in the plot
    nUpdates (1,1) {mustBeInteger, mustBePositive} = 10 % Number of times data is extracted for plots during each pulse. A minimum of 1 update is performed in each phase (2 for each pulse consisting of an illumination phase and a diffusion phase)
    mediaPropRecalcPeriod (1,1) {mustBeInteger, mustBePositive} = 1 % Every N updates, the media properties will be recalculated (including, if needed, re-running MC and FMC steps)

    slicePositions (1,3) {mustBeInRange(slicePositions,0,1)} = [.5 1 1] % Relative slice positions [x y z] for the 3D plots on a scale from 0 to 1
    tempSensorPositions (:,3) {mustBeFinite} % Each row is a temperature sensor's absolute [x y z] coordinates. Leave the matrix empty ([]) to disable temperature sensors.

    %% Calculated properties
    mediaProperties = NaN
    mediaProperties_funcHandles = NaN

    M = NaN

    Tdependent = NaN
    FDdependent = NaN

    T  = NaN                                % Final temperature
    Omega single = NaN
    maxMediaTemps = NaN
    sensorsTimeVector = NaN
    sensorTemps = NaN
    movieFrames = NaN
  end

  methods
    function obj = heatSimulation()
      %HEATSIMULATION Construct an instance of this class
    end
  end
end

function mustBeAllFiniteOrAllNaN(x)
x = x(:);
if all(isnan(x)) || all(isfinite(x))
  % Valid input
else
  error('Value must be all NaNs or all finite values.');
end
end

function mustBeFiniteOrNaNScalar(x)
x = x(:);
if isscalar(x)
  if isnan(x) || isfinite(x)
    % Valid input
  else
    error('Value must be either a scalar NaN or a finite scalar or array.');
  end
else
  if all(isfinite(x))
    % Valid input
  else
    error('Value must be either a scalar NaN or a finite scalar or array.');
  end
end
end

function mustBeInRange(x,a,b)
if any(x(:) < a) || any(x(:) > b)
  error('Error: Values must be in range %f to %f.',a,b);
end
end