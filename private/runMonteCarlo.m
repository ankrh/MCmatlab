function MCoutput = runMonteCarlo(MCinput)
%   Created 2018 by Dominik Marti and Anders K. Hansen, DTU Fotonik
%
%   Prepares the illumination beam and runs the Monte Carlo simulation.
%   After finishing, calls plotMCmatlab for the display of the result.
%
%	Pay attention to the sections with headers that say "USER SPECIFIED:"
%	In those sections, you must fill in the parameters relevant for your simulation.
%
%   Input
%       name
%           the basename of the file saved by defineGeometry.m
%
%   Output
%       ./Data/[name]_MCoutput.mat
%           file containing the 3D fluence rate distribution
%
%   Requires
%       deleteDataFiles.m
%       MCmatlab.mex (architecture specific)
%       plotMCmatlab.m
%

%% Check to ensure that the light collector is not inside the cuboid and set res to 1 if using fiber
if isfield(MCinput,'LightCollector')
	MCinput.useLightCollector = true;
    if isfinite(MCinput.LightCollector.f)
        xLCC = MCinput.LightCollector.xFPC - MCinput.LightCollector.f*sin(MCinput.LightCollector.theta)*cos(MCinput.LightCollector.phi); % x position of Light Collector Center
        yLCC = MCinput.LightCollector.yFPC - MCinput.LightCollector.f*sin(MCinput.LightCollector.theta)*sin(MCinput.LightCollector.phi); % y position
        zLCC = MCinput.LightCollector.zFPC - MCinput.LightCollector.f*cos(MCinput.LightCollector.theta);             % z position
    else
        xLCC = MCinput.LightCollector.xFPC;
        yLCC = MCinput.LightCollector.yFPC;
        zLCC = MCinput.LightCollector.zFPC;
        MCinput.LightCollector.res = 1;
    end

    if (abs(xLCC)                           < MCinput.G.nx*MCinput.G.dx/2 && ...
        abs(yLCC)                           < MCinput.G.ny*MCinput.G.dy/2 && ...
        abs(zLCC - MCinput.G.nz*MCinput.G.dz/2) < MCinput.G.nz*MCinput.G.dz/2)
        error('Error: Light collector center (%.4f,%.4f,%.4f) is inside cuboid',xLCC,yLCC,zLCC);
    end

	if ~isfield(MCinput.LightCollector,'tStart')
		MCinput.LightCollector.tStart = 0;
		MCinput.LightCollector.tEnd = 0;
		MCinput.LightCollector.nTimeBins = 0;
	end
else
	MCinput.useLightCollector = false;
	MCinput.LightCollector.xFPC = 0;
	MCinput.LightCollector.yFPC = 0;
	MCinput.LightCollector.zFPC = 0;
	MCinput.LightCollector.theta = 0;
	MCinput.LightCollector.phi = 0;
	MCinput.LightCollector.f = 0;
	MCinput.LightCollector.diam = 0;
	MCinput.LightCollector.FieldSize = 0;
	MCinput.LightCollector.NA = 0;
	MCinput.LightCollector.res = 0;
	MCinput.LightCollector.tStart = 0;
	MCinput.LightCollector.tEnd = 0;
	MCinput.LightCollector.nTimeBins = 0;
end

%% Call Monte Carlo C script (MEX file) to get fluence rate (intensity) distribution
MCinput.G.M = MCinput.G.M-1; % Convert to C-style indexing
MCoutput = MCmatlab(MCinput);
clear MCmatlab; % Unload MCmatlab MEX file so it can be modified externally again

if(~MCinput.silentMode)
	MCinput.G.M = MCinput.G.M+1; % Convert back to MATLAB-style indexing
	plotMCmatlab(MCinput,MCoutput);
end
end
