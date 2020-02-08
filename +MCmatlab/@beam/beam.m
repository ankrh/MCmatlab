classdef beam
    %BEAM Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        beamType                                % 0: Pencil beam, 1: Isotropically emitting point source, 2: Infinite plane wave, 3: Laguerre-Gaussian LG01 beam, 4: Radial-factorizable beam (e.g., a Gaussian beam), 5: X/Y factorizable beam (e.g., a rectangular LED emitter)
        
        xFocus = 0
        yFocus = 0
        zFocus = 0
        
        theta  = 0                              % [rad] Polar angle of beam center axis
        phi  = 0                                % [rad] Azimuthal angle of beam center axis
        psi  = 0                                % [rad] Axial rotation angle of beam, relevant only for XY distributed beams

        NF struct
        FF struct
    end
    
    properties (Hidden)
        sourceDistribution
    end
    
    methods
        function obj = beam()
            %BEAM Construct an instance of this class
            %   Detailed explanation goes here
            
        end
    end
end

