classdef lightCollector
  %LIGHTCOLLECTOR This class contains all properties related to a
  %light collector of an MCmatlab.model.
  %   This class defines the properties of a lightcollector to be used in
  %   monteCarloSimulation or fluorescenceMonteCarloSimulation.

  properties
    %% Input properties
    x (1,1) double {mustBeFinite} = 0 % [cm] x position of either the center of the objective lens focal plane or the fiber tip
    y (1,1) double {mustBeFinite} = 0 % [cm] y position
    z (1,1) double {mustBeFinite} = 0 % [cm] z position

    theta (1,1) double {mustBeFinite} = 0 % [rad] Polar angle of direction the light collector is facing
    phi (1,1) double {mustBeFinite} = pi/2 % [rad] Azimuthal angle of direction the light collector is facing

    f (1,1) double {mustBePositive} = Inf % [cm] Focal length of the objective lens (if light collector is a fiber, set this to Inf).
    diam (1,1) double {mustBePositive} = 1 % [cm] Diameter of the light collector aperture. For an ideal thin lens, this is 2*f*tan(asin(NA)).
    fieldSize (1,1) double {mustBePositive} = 1 % [cm] Field Size of the imaging system (diameter of area in object plane that gets imaged). Only used for finite f.
    NA (1,1) double {mustBePositive} = 0.22 % [-] Fiber NA. Only used for infinite f.
    res (1,1) double {mustBeInteger, mustBeNonnegative} = 0 % X and Y resolution of light collector in pixels, only used for finite f

    tStart (1,1) double {mustBeFinite} = 0 % [s] Start of the detection time interval
    tEnd (1,1) double {mustBeFinite} = 1e-12 % [s] End of the detection time interval
    nTimeBins (1,1) double {mustBeInteger, mustBeNonnegative} = 0 % Number of bins between tStart and tEnd. If zero, the measurement is not time-resolved.

    %% Calculated properties
    image = NaN
    X = NaN
    Y = NaN
    t = NaN
  end

  methods
    function obj = lightCollector()
      %LIGHTCOLLECTOR Construct an instance of this class

    end
  end
end

