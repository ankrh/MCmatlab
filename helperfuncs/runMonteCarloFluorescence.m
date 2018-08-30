function FMCoutput = runMonteCarloFluorescence(FMCinput)
%   Created 2018 by Dominik Marti and Anders K. Hansen, DTU Fotonik
%
%   Script for simulating distribution and magnitude of fluorescence
%   based on the output of runMonteCarlo.m
%
%   Prepares and runs the Monte Carlo simulation.
%   After finishing, calls plotMCmatlab for the display of the result.
%
%	Pay attention to the sections with headers that say "USER SPECIFIED:"
%	In those sections, you must fill in the parameters relevant for your simulation.
%
%   Input
%       name
%           the basename of the files saved by defineGeometry.m and runMonteCarlo.m
%
%   Output
%       ./Data/[name]_MCoutput_fluorescence.mat
%           file containing the 3D fluorescence fluence rate distribution
%
%   Requires
%       deleteDataFiles.m
%       MCmatlab.mex (architecture specific)
%       plotMCmatlab.m
%

%% Calculate 3D fluorescence source distribution, including saturation
mua_vec = [FMCinput.G.mediaProperties.mua]; % The media's excitation absorption coefficients
Y_vec = [FMCinput.G.mediaProperties.Y]; % The media's fluorescence power yields
sat_vec = [FMCinput.G.mediaProperties.sat]; % The media's fluorescence saturation fluence rates (intensity)
FMCinput.Beam.sourceDistribution = Y_vec(FMCinput.G.M).*mua_vec(FMCinput.G.M)*FMCinput.Beam.P_excitation.*FMCinput.MCoutput.F./(1 + FMCinput.Beam.P_excitation*FMCinput.MCoutput.F./sat_vec(FMCinput.G.M)); % [W/cm^3]
if(max(FMCinput.Beam.sourceDistribution(:)) == 0); error('Error: No fluorescence emitters'); end

%% Check to ensure that the light collector is not inside the cuboid and set res_LC to 1 if using fiber
if isfield(FMCinput,'LightCollector')
	FMCinput.useLightCollector = true;
    if isfinite(FMCinput.LightCollector.f)
        xLCC = FMCinput.LightCollector.x - FMCinput.LightCollector.f*sin(FMCinput.LightCollector.theta)*cos(FMCinput.LightCollector.phi); % x position of Light Collector Center
        yLCC = FMCinput.LightCollector.y - FMCinput.LightCollector.f*sin(FMCinput.LightCollector.theta)*sin(FMCinput.LightCollector.phi); % y position
        zLCC = FMCinput.LightCollector.z - FMCinput.LightCollector.f*cos(FMCinput.LightCollector.theta);             % z position
    else
        xLCC = FMCinput.LightCollector.x;
        yLCC = FMCinput.LightCollector.y;
        zLCC = FMCinput.LightCollector.z;
        FMCinput.LightCollector.res = 1;
    end

    if (abs(xLCC)                           < FMCinput.G.nx*FMCinput.G.dx/2 && ...
        abs(yLCC)                           < FMCinput.G.ny*FMCinput.G.dy/2 && ...
        abs(zLCC - FMCinput.G.nz*FMCinput.G.dz/2) < FMCinput.G.nz*FMCinput.G.dz/2)
        error('Error: Light collector center (%.4f,%.4f,%.4f) is inside cuboid',xLCC,yLCC,zLCC);
    end

	if ~isfield(FMCinput.LightCollector,'tStart')
		FMCinput.LightCollector.tStart = 0;
		FMCinput.LightCollector.tEnd = 0;
		FMCinput.LightCollector.nTimeBins = 0;
	end
else
	FMCinput.useLightCollector = false;
	FMCinput.LightCollector.x = 0;
	FMCinput.LightCollector.y = 0;
	FMCinput.LightCollector.z = 0;
	FMCinput.LightCollector.theta = 0;
	FMCinput.LightCollector.phi = 0;
	FMCinput.LightCollector.f = 0;
	FMCinput.LightCollector.diam = 0;
	FMCinput.LightCollector.FieldSize = 0;
	FMCinput.LightCollector.NA = 0;
	FMCinput.LightCollector.res = 0;
	FMCinput.LightCollector.tStart = 0;
	FMCinput.LightCollector.tEnd = 0;
	FMCinput.LightCollector.nTimeBins = 0;
end

FMCinput.G.M = FMCinput.G.M - 1; % Convert to C-style indexing
FMCoutput = MCmatlab(FMCinput); % FMCoutput.F is an absolute fluence rate (intensity) quantity, unlike the non-fluorescence MCoutput.F which are actually fluence rates normalized to the incident power
clear MCmatlab; % Unload MCmatlab MEX file so it can be modified externally again

if(~FMCinput.silentMode)
    FMCinput.G.M = FMCinput.G.M + 1; % Convert back to MATLAB-style indexing
    plotMCmatlabFluorescence(FMCinput,FMCoutput);
end
end
