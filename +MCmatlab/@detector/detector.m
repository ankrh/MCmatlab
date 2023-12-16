classdef detector
  %   This class defines the properties of a detector to be used in
  %   monteCarloSimulation or fluorescenceMonteCarloSimulation.

  properties
    %% Input properties
    % Position coordinates of the center of the detector face:
    x (1,1) double {mustBeFinite} = 0 % [cm] x position
    y (1,1) double {mustBeFinite} = 0 % [cm] y position
    z (1,1) double {mustBeFinite} = 0 % [cm] z position

    theta (1,1) double {mustBeFinite} = 0 % [rad] Polar angle of direction the detector is facing
    phi (1,1) double {mustBeFinite} = pi/2 % [rad] Azimuthal angle of direction the detector is facing
    psi (1,1) double {mustBeFinite} = 0 % [rad] Axial rotation angle of the detector, determining the direction of the X and Y axes. Relevant only if Xsize ~= Ysize or if you want an image to have a particular orientation

    Xsize (1,1) double {mustBePositive} = 1 % Extent of the detector in the X direction (not the same as the x direction)
    Ysize (1,1) double {mustBePositive} = 1 % Extent of the detector in the Y direction (not the same as the y direction)

    shape (1,1) MCmatlab.shape = 'Rectangle' % Can be either ellipse or rectangle.
    
    res (1,1) double {mustBeInteger, mustBeNonnegative} = 1 % X and Y resolution of detector in pixels

    tStart (1,1) double {mustBeFinite} = 0 % [s] Start of the detection time interval
    tEnd (1,1) double {mustBeFinite} = 1e-12 % [s] End of the detection time interval
    nTimeBins (1,1) double {mustBeInteger, mustBeNonnegative} = 0 % Number of bins between tStart and tEnd. If zero, the measurement is not time-resolved.

    %% Calculated properties
    irradiance = NaN % (X,Y,lambda) array
    power = NaN % Scalar, integral of irradiance
    X = NaN
    Y = NaN
    t = NaN
  end

  properties (Hidden, Dependent)
    image
  end

  methods
    function x   = get.image(obj  ); x = obj.irradiance; end
    function obj = set.image(obj,x);     obj.irradiance = x; end  
  end
end

