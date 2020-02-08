classdef lightCollector
    %LIGHTCOLLECTOR Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        x                                       % [cm] x position of either the center of the objective lens focal plane or the fiber tip
        y                                       % [cm] y position
        z                                       % [cm] z position
        
        theta  = 0                              % [rad] Polar angle of direction the light collector is facing
        phi  = pi/2                             % [rad] Azimuthal angle of direction the light collector is facing
        
        f                                       % [cm] Focal length of the objective lens (if light collector is a fiber, set this to Inf).
        diam                                    % [cm] Diameter of the light collector aperture. For an ideal thin lens, this is 2*f*tan(asin(NA)).
        fieldSize                               % [cm] Field Size of the imaging system (diameter of area in object plane that gets imaged). Only used for finite f.
        NA                                      % [-] Fiber NA. Only used for infinite f.
        res                                     % X and Y resolution of light collector in pixels, only used for finite f

        tStart                                  % [s] Start of the detection time interval
        tEnd                                    % [s] End of the detection time interval
        nTimeBins  = 0                          % Number of bins between tStart and tEnd. If zero, the measurement is not time-resolved.
    end
    
    properties (Hidden)
        image
        X
        Y
    end
    
    methods
        function obj = lightCollector()
            %LIGHTCOLLECTOR Construct an instance of this class
            %   Detailed explanation goes here
            
        end
    end
end

