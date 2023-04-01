classdef fluorescenceMonteCarloSimulation
  %FLUORESCENCEMONTECARLOSIMULATION This class includes all properties and methods
  %related to the fluorescence Monte Carlo simulation in an MCmatlab.model.

  properties
    %% Input properties
    silentMode (1,1) logical = false % Disables command window text and progress indication
    useAllCPUs (1,1) logical = false % If false, MCmatlab will leave one processor unused. Useful for doing other work on the PC while simulations are running.
    useGPU (1,1) logical = false % Use CUDA acceleration for NVIDIA GPUs
    GPUdevice (1,1) int32 = 0 % GPU device index to run on (default: 0, the first one)
    simulationTimeRequested (1,1) double {mustBePositive} = 0.1 % [min] Time duration of the simulation
    nPhotonsRequested (1,1) double {mustBeFinitePositiveIntegerOrNaN} = NaN % # of photons to launch
    requestCollectedPhotons (1,1) logical = false % If true, the photon # in nPhotonsRequested is interpreted as collected photons rather than launched photons
    calcNormalizedFluenceRate (1,1) logical = true % If true, the 3D normalized fluence rate output array will be calculated. Set to false if you have a light collector and you're only interested in the image output.
    nExamplePaths (1,1) double {mustBeInteger, mustBeNonnegative} = 0 % This number of photons will have their paths stored and shown after completion, for illustrative purposes
    farFieldRes (1,1) double {mustBeInteger, mustBeNonnegative} = 0 % If nonzero, photons that "escape" will have their energies tracked in a 2D angle distribution (theta,phi) array with theta and phi resolutions equal to this number. An "escaping" photon is one that hits the top cuboid boundary (if boundaryType == 2) or any cuboid boundary (if boundaryType == 1) where the medium has refractive index 1.

    matchedInterfaces (1,1) logical = true % If true, assumes all refractive indices are 1. If false, uses the refractive indices defined in getMediaProperties
    smoothingLengthScale (1,1) double {mustBePositive} = 0.1 % Length scale over which smoothing of the Sobel interface gradients should be performed
    boundaryType (1,1) double {mustBeInteger, mustBeInRange(boundaryType,0,3)} = 1 % 0: No escaping boundaries, 1: All cuboid boundaries are escaping, 2: Top cuboid boundary only is escaping, 3: Top and bottom boundaries are escaping, while the side boundaries are cyclic
    wavelength (1,:) double {mustBeFinitePositiveOrNaN} = NaN % [nm] Excitation wavelength, used for determination of optical properties for excitation light

    useLightCollector (1,1) logical = false
    lightCollector (1,1) MCmatlab.lightCollector

    depositionCriteria (1,1) MCmatlab.depositionCriteria

    %% Calculated properties
    simulationTime = NaN;
    nPhotons = NaN;
    nPhotonsCollected = NaN
    nThreads = NaN;

    mediaProperties = NaN; % Wavelength- and splitting-dependent
    CDFs = NaN % Cumulative distribution functions for custom phase functions (not Henyey Greenstein)
    M = NaN; % Splitting-dependent
    interfaceNormals single = NaN

    examplePaths = NaN;

    normalizedFluenceRate single = 0
    normalizedFluenceRate_detected = NaN

    farField = NaN;
    farFieldTheta = NaN;
    farFieldPhi = NaN;

    sourceDistribution = NaN;

    normalizedIrradiance_xpos = NaN % Normalized irradiance on the boundary in the positive x direction
    normalizedIrradiance_xneg = NaN
    normalizedIrradiance_ypos = NaN
    normalizedIrradiance_yneg = NaN
    normalizedIrradiance_zpos = NaN
    normalizedIrradiance_zneg = NaN
  end

  properties (Dependent)
    normalizedAbsorption
    normalizedAbsorption_detected
  end

  properties (Hidden)
    calcNFR
    calcNFRdet
    LC
    NFR
    NFRdet
    NA
    NAdet
    NI_xpos
    NI_xneg
    NI_ypos
    NI_yneg
    NI_zpos
    NI_zneg
  end

  methods
    function NA = get.normalizedAbsorption(obj)
      if numel(obj.NFR) > 1
        NA = NaN(size(obj.NFR));
        for iWavelength = 1:numel(obj.wavelength)
          mua_vec = [obj.mediaProperties.mua(:,iWavelength)];
          NA(:,:,:,iWavelength) = mua_vec(obj.M).*obj.NFR(:,:,:,iWavelength);
        end
      else
        NA = 0;
      end
    end
    function NAdet = get.normalizedAbsorption_detected(obj)
      if numel(obj.NFRdet) > 1
        NAdet = NaN(size(obj.NFRdet));
        for iWavelength = 1:size(obj.mediaProperties,2)
          mua_vec = [obj.mediaProperties.mua(:,iWavelength)];
          NAdet(:,:,:,iWavelength) = mua_vec(obj.M).*obj.NFRdet(:,:,:,iWavelength);
        end
      else
        NAdet = 0;
      end
    end
    function x   = get.calcNFR(obj  ); x = obj.calcNormalizedFluenceRate; end
    function obj = set.calcNFR(obj,x);     obj.calcNormalizedFluenceRate = x; end %#ok<MCSUP> 
    function x   = get.calcNFRdet(obj  ); x = obj.calcNormalizedFluenceRate_detected; end
    function obj = set.calcNFRdet(obj,x);     obj.calcNormalizedFluenceRate_detected = x; end %#ok<MCSUP> 
    function x   = get.LC(obj  ); x = obj.lightCollector; end
    function obj = set.LC(obj,x);     obj.lightCollector = x; end %#ok<MCSUP> 
    function x   = get.NFR(obj  ); x = obj.normalizedFluenceRate; end
    function obj = set.NFR(obj,x);     obj.normalizedFluenceRate = x; end %#ok<MCSUP> 
    function x   = get.NFRdet(obj  ); x = obj.normalizedFluenceRate_detected; end
    function obj = set.NFRdet(obj,x);     obj.normalizedFluenceRate_detected = x; end %#ok<MCSUP> 
    function x   = get.NA(obj  ); x = obj.normalizedAbsorption; end
    function obj = set.NA(obj,x);     obj.normalizedAbsorption = x; end %#ok<MCSUP> 
    function x   = get.NAdet(obj  ); x = obj.normalizedAbsorption_detected; end
    function obj = set.NAdet(obj,x);     obj.normalizedAbsorption_detected = x; end %#ok<MCSUP> 
    function x   = get.NI_xpos(obj  ); x = obj.normalizedIrradiance_xpos; end
    function obj = set.NI_xpos(obj,x);     obj.normalizedIrradiance_xpos = x; end %#ok<MCSUP> 
    function x   = get.NI_xneg(obj  ); x = obj.normalizedIrradiance_xneg; end
    function obj = set.NI_xneg(obj,x);     obj.normalizedIrradiance_xneg = x; end %#ok<MCSUP> 
    function x   = get.NI_ypos(obj  ); x = obj.normalizedIrradiance_ypos; end
    function obj = set.NI_ypos(obj,x);     obj.normalizedIrradiance_ypos = x; end %#ok<MCSUP> 
    function x   = get.NI_yneg(obj  ); x = obj.normalizedIrradiance_yneg; end
    function obj = set.NI_yneg(obj,x);     obj.normalizedIrradiance_yneg = x; end %#ok<MCSUP> 
    function x   = get.NI_zpos(obj  ); x = obj.normalizedIrradiance_zpos; end
    function obj = set.NI_zpos(obj,x);     obj.normalizedIrradiance_zpos = x; end %#ok<MCSUP> 
    function x   = get.NI_zneg(obj  ); x = obj.normalizedIrradiance_zneg; end
    function obj = set.NI_zneg(obj,x);     obj.normalizedIrradiance_zneg = x; end %#ok<MCSUP> 
  end
end

function mustBeFinitePositiveIntegerOrNaN(x)
x = x(:);
if all(isnan(x) | (isfinite(x) & rem(x,1) == 0 & x > 0))
  % Input valid
else
  error('Value must be a finite positive integer or NaN.');
end
end

function mustBeFinitePositiveOrNaN(x)
x = x(:);
if all(isnan(x) | (isfinite(x) & x > 0))
  % Valid input
else
  error('Value must be finite positive or NaN.');
end
end

function mustBeInRange(x,a,b)
if any(x(:) < a) || any(x(:) > b)
  error('Error: Values must be in range %f to %f.',a,b);
end
end

function mustBeScalarOrEmpty(x)
if numel(x) > 1
  error('Error: Value must be scalar or empty.');
end
end