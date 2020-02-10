classdef beamField
    %BEAMFIELD Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        radialDistr % Radial near or far field distribution - 0: Top-hat, 1: Gaussian, Array: Custom. Doesn't need to be normalized.
        radialWidth % [cm] Radial near field 1/e^2 radius if top-hat or Gaussian or half-width of the full distribution if custom
        % [rad] Radial far field 1/e^2 half-angle if top-hat or Gaussian or half-angle of the full distribution if custom. For a diffraction limited Gaussian beam, this should be set to model.MC.wavelength*1e-9/(pi*model.MC.beam.NF.radialWidth*1e-2))
        XDistr % X near or far field distribution - 0: Top-hat, 1: Gaussian, Array: Custom. Doesn't need to be normalized.
        XWidth % [cm] X near field 1/e^2 radius if top-hat or Gaussian or half-width of the full distribution if custom
        % [rad] X far field 1/e^2 half-angle if top-hat or Gaussian or half-angle of the full distribution if custom
        YDistr % Y near or far field distribution - 0: Top-hat, 1: Gaussian, Array: Custom. Doesn't need to be normalized.
        YWidth % [cm] Y near field 1/e^2 radius if top-hat or Gaussian or half-width of the full distribution if custom
        % [rad] Y far field 1/e^2 half-angle if top-hat or Gaussian or half-angle of the full distribution if custom

    end
    
    methods
        function obj = beamField()
            %BEAMFIELD Construct an instance of this class
            %   Detailed explanation goes here
            
        end
    end
end

