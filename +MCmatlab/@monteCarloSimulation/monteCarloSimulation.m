classdef monteCarloSimulation
    %monteCarloSimulation This class includes all properties and methods
    %related to the Monte Carlo simulation in an MCmatlab.model.
    
    properties
        silentMode logical = false              % Disables command window text and progress indication
        useAllCPUs logical = false              % If false, MCmatlab will leave one processor unused. Useful for doing other work on the PC while simulations are running.
        useGPU logical = false                  % Use CUDA acceleration for NVIDIA GPUs
        simulationTimeRequested = 0.1           % [min] Time duration of the simulation
        nPhotonsRequested = NaN                 % # of photons to launch
        calcNFR logical = true                  % If true, the 3D fluence rate output array NFR will be calculated. Set to false if you have a light collector and you're only interested in the image output.
        calcNFRdet logical = false              % If true, the 3D fluence rate output array NFRdet will be calculated. Only photons that end up on the light collector are counted in NFRdet.
        nExamplePaths = 0                       % This number of photons will have their paths stored and shown after completion, for illustrative purposes
        farFieldRes = 0                         % If nonzero, photons that "escape" will have their energies tracked in a 2D angle distribution (theta,phi) array with theta and phi resolutions equal to this number. An "escaping" photon is one that hits the top cuboid boundary (if boundaryType == 2) or any cuboid boundary (if boundaryType == 1) where the medium has refractive index 1.
        matchedInterfaces logical = true        % If true, assumes all refractive indices are 1. If false, uses the refractive indices defined in getMediaProperties
        interpolateNormals logical = true       % If true, will perform linear 3D interpolation to make a more precise estimate of the local surface gradient at interfaces
        smoothingLengthScale = 0                % Length scale over which smoothing of the Sobel interface gradients should be performed
        boundaryType = 1                        % 0: No escaping boundaries, 1: All cuboid boundaries are escaping, 2: Top cuboid boundary only is escaping
        wavelength = NaN                        % [nm] Excitation wavelength, used for determination of optical properties for excitation light
        beam MCmatlab.beam
        P = NaN                                 % [W] Incident pulse peak power (in case of infinite plane waves, only the power incident upon the cuboid's top surface)
        FRinitial = NaN                         % [W/cm^2] Initial guess for the intensity distribution, to be used for fluence rate dependent simulations
        FR = NaN
        FRdepIterations = 20
        useLightCollector logical = false
        LC MCmatlab.lightCollector
        
        simulationTime = NaN
        nPhotons = NaN
        nThreads = NaN

        mediaProperties_funcHandles = NaN % Wavelength-dependent
        mediaProperties = NaN % Wavelength- and splitting-dependent
        FRdependent = NaN
        FDdependent = NaN
        Tdependent = NaN
        M = NaN % Splitting-dependent
        interfaceNormals = NaN

        examplePaths = NaN

        NFR = NaN % Normalized Fluence Rate
        NFRdet = NaN

        farField = NaN
        farFieldTheta = NaN
        farFieldPhi = NaN

        NI_xpos = NaN % Normalized irradiance on the boundary in the positive x direction
        NI_xneg = NaN
        NI_ypos = NaN
        NI_yneg = NaN
        NI_zpos = NaN
        NI_zneg = NaN
    end

    
    methods
        function obj = monteCarloSimulation()
            %monteCarloSimulation Construct an instance of this class
            
            obj.beam = MCmatlab.beam;
            obj.LC = MCmatlab.lightCollector;
        end

    end
end

