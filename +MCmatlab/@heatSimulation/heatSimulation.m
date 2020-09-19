classdef heatSimulation
    %HEATSIMULATION Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        silentMode logical = false              % Disables command window text and progress indication
        useAllCPUs logical = false              % If false, MCmatlab will leave one processor unused. Useful for doing other work on the PC while simulations are running.
        useGPU logical = false                  % Use CUDA acceleration for NVIDIA GPUs
        makeMovie logical = false               % Requires silentMode = false.
        deferMovieWrite logical = false
        largeTimeSteps logical = false          % If true, calculations will be faster, but some voxel temperatures may be slightly less precise. Test for yourself whether this precision is acceptable for your application.
        
        heatBoundaryType = 0                    % 0: Insulating boundaries, 1: Constant-temperature boundaries (heat-sinked)
        durationOn = NaN                        % [s] Pulse on-duration
        durationOff = NaN                       % [s] Pulse off-duration
        durationEnd = NaN                       % [s] Non-illuminated relaxation time to add to the end of the simulation to let temperature diffuse after the pulse train
        Tinitial = NaN                          % [deg C] Initial temperature, can be scalar or 3D array
        
        nPulses = 1                             % Number of consecutive pulses, each with an illumination phase and a diffusion phase. If simulating only illumination or only diffusion, use n_pulses = 1.
        
        plotTempLimits = NaN                    % [deg C] Expected range of temperatures, used only for setting the color scale in the plot
        nUpdates = 10                           % Number of times data is extracted for plots during each pulse. A minimum of 1 update is performed in each phase (2 for each pulse consisting of an illumination phase and a diffusion phase)
        mediaPropRecalcPeriod = 1               % Every N updates, the media properties will be recalculated (including, if needed, re-running MC and FMC steps)
        
        slicePositions = [.5 1 1]               % Relative slice positions [x y z] for the 3D plots on a scale from 0 to 1
        tempSensorPositions                     % Each row is a temperature sensor's absolute [x y z] coordinates. Leave the matrix empty ([]) to disable temperature sensors. 

        mediaProperties = NaN
        mediaProperties_funcHandles = NaN
        
        M = NaN
        
        Tdependent = NaN
        FDdependent = NaN

        T  = NaN                                % Final temperature
        Omega single = NaN
        maxMediaTemps = NaN
        sensorsTimeVector = NaN
        sensorTemps = NaN
        movieFrames = NaN

    end
    
    methods
        function obj = heatSimulation()
            %HEATSIMULATION Construct an instance of this class
            %   Detailed explanation goes here
            
        end

    end
end

