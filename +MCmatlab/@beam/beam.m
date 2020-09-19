classdef beam
    %BEAM Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        beamType = NaN                          % 0: Pencil beam, 1: Isotropically emitting point source, 2: Infinite plane wave, 3: Laguerre-Gaussian LG01 beam, 4: Radial-factorizable beam (e.g., a Gaussian beam), 5: X/Y factorizable beam (e.g., a rectangular LED emitter)
        
        xFocus = 0
        yFocus = 0
        zFocus = 0
        
        theta  = 0                              % [rad] Polar angle of beam center axis
        phi  = 0                                % [rad] Azimuthal angle of beam center axis
        psi  = 0                                % [rad] Axial rotation angle of beam, relevant only for XY distributed beams

        NF MCmatlab.beamField
        FF MCmatlab.beamField

        sourceDistribution = NaN
    end
    
    methods
        function obj = beam()
            %BEAM Construct an instance of this class
            %   Detailed explanation goes here
            
            obj.NF = MCmatlab.beamField;
            obj.FF = MCmatlab.beamField;
        end
    end
end

