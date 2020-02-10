classdef fluorescenceMonteCarloSimulation
    %FLUORESCENCEMONTECARLOSIMULATION Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        useGPU = false;

        simulationTimeRequested = 0.1;
        nPhotonsRequested = NaN;
        silentMode = false;
        useAllCPUs = false;
        calcNFR = true;
        calcNFRdet = false;
        nExamplePaths = 0;
        farFieldRes = 0;

        matchedInterfaces = true;
        boundaryType = 1;
        wavelength = NaN;
    end
    
    properties (Hidden)
        %% Fluorescence Monte Carlo parameters that are calculated
        simulationTime = NaN;
        nPhotons = NaN;
        nThreads = NaN;

        mediaProperties_funcHandles = NaN; % Wavelength-dependent
        mediaProperties = NaN; % Wavelength- and splitting-dependent
        FRdependent = NaN;
        FDdependent = NaN;
        Tdependent = NaN;
        M = NaN; % Splitting-dependent
        RI = NaN;

        examplePaths = NaN;

        NFR = NaN; % Normalized Fluence Rate
        NFRdet = NaN;

        farField = NaN;
        farFieldTheta = NaN;
        farFieldPhi = NaN;

        sourceDistribution = NaN;

        NI_xpos = NaN; % Normalized irradiance on the boundary in the positive x direction
        NI_xneg = NaN;
        NI_ypos = NaN;
        NI_yneg = NaN;
        NI_zpos = NaN;
        NI_zneg = NaN;
    end
    
    methods
        function obj = fluorescenceMonteCarloSimulation()
            %FLUORESCENCEMONTECARLOSIMULATION Construct an instance of this class
            %   Detailed explanation goes here
            
        end
        
    end
end

