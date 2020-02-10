classdef geometry
    %GEOMETRY Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        silentMode logical = false              % Disables command window text and progress indication
        nx = NaN                                % Number of bins in the x direction
        ny = NaN                                % Number of bins in the y direction
        nz = NaN                                % Number of bins in the z direction
        Lx = NaN                                % [cm] x size of simulation cuboid
        Ly = NaN                                % [cm] y size of simulation cuboid
        Lz = NaN                                % [cm] z size of simulation cuboid
        mediaPropertiesFunc function_handle     % Media properties defined as a function at the end of the model file
        mediaPropParams cell = {}               % Cell array containing any additional parameters to be passed to the getMediaProperties function
        geomFunc function_handle                % Function to use for defining the distribution of media in the cuboid. Defined at the end of the model file.
        geomFuncParams cell = {}                % Cell array containing any additional parameters to pass into the geometry function, such as media depths, inhomogeneity positions, radii etc.
    end
    
    properties (Hidden)
        dx = NaN
        dy = NaN
        dz = NaN
        x = NaN
        y = NaN
        z = NaN
        M_raw uint8 = NaN
    end
    
    methods
        function obj = geometry()
            %GEOMETRY Construct an instance of this class
            %   Detailed explanation goes here
            
        end
                
    end
end

