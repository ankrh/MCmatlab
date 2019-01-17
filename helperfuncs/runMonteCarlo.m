function MCoutput = runMonteCarlo(MCinput)
%   Requires
%       MCmatlab.mex (architecture specific)
%
%	See also defineGeometry, plotMCmatlab, runMonteCarloFluorescence, simulateHeatDistribution

%%%%%
%   Copyright 2017, 2018 by Dominik Marti and Anders K. Hansen, DTU Fotonik
%
%   This file is part of MCmatlab.
%
%   MCmatlab is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
%
%   MCmatlab is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with MCmatlab.  If not, see <https://www.gnu.org/licenses/>.
%%%%%

if ~isfield(MCinput,'silentMode')
	MCinput.silentMode = false;
end
if ~isfield(MCinput,'useAllCPUs')
	MCinput.useAllCPUs = false;
end
if ~isfield(MCinput,'calcF')
    MCinput.calcF = true;
end
if ~isfield(MCinput,'calcFdet')
    MCinput.calcFdet = false;
end

if isfield(MCinput,'LightCollector')
    %% Check to ensure that the light collector is not inside the cuboid and set res to 1 if using fiber
	MCinput.useLightCollector = true;
    if isfinite(MCinput.LightCollector.f)
        xLCC = MCinput.LightCollector.x - MCinput.LightCollector.f*sin(MCinput.LightCollector.theta)*cos(MCinput.LightCollector.phi); % x position of Light Collector Center
        yLCC = MCinput.LightCollector.y - MCinput.LightCollector.f*sin(MCinput.LightCollector.theta)*sin(MCinput.LightCollector.phi); % y position
        zLCC = MCinput.LightCollector.z - MCinput.LightCollector.f*cos(MCinput.LightCollector.theta);             % z position
    else
        xLCC = MCinput.LightCollector.x;
        yLCC = MCinput.LightCollector.y;
        zLCC = MCinput.LightCollector.z;
        MCinput.LightCollector.res = 1;
    end

    if (abs(xLCC)                           < MCinput.G.nx*MCinput.G.dx/2 && ...
        abs(yLCC)                           < MCinput.G.ny*MCinput.G.dy/2 && ...
        abs(zLCC - MCinput.G.nz*MCinput.G.dz/2) < MCinput.G.nz*MCinput.G.dz/2)
        error('Error: Light collector center (%.4f,%.4f,%.4f) is inside cuboid',xLCC,yLCC,zLCC);
    end

    %% If no time tagging start value was defined, assume no time tagging is to be performed
	if ~isfield(MCinput.LightCollector,'tStart')
		MCinput.LightCollector.tStart = 0;
		MCinput.LightCollector.tEnd = 0;
		MCinput.LightCollector.nTimeBins = 0;
	end
else
    %% Assume no light collector is present
	MCinput.useLightCollector = false;
	MCinput.LightCollector.x = 0;
	MCinput.LightCollector.y = 0;
	MCinput.LightCollector.z = 0;
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

if ~MCinput.calcF && ~MCinput.useLightCollector
    error('Error: calcF is false, but no light collector is defined');
end
if MCinput.calcFdet && ~MCinput.useLightCollector
    error('Error: calcFdet is true, but no light collector is defined');
end

%% Call Monte Carlo C script (MEX file) to get fluence rate (intensity) distribution
MCinput.G.M = MCinput.G.M-1; % Convert to C-style indexing
MCoutput = MCmatlab(MCinput);

% Add positions of the centers of the pixels in the light collector image
if MCinput.useLightCollector && MCinput.LightCollector.res > 1
	MCoutput.X = linspace(MCinput.LightCollector.FieldSize*(1/MCinput.LightCollector.res-1),MCinput.LightCollector.FieldSize*(1-1/MCinput.LightCollector.res),MCinput.LightCollector.res)/2;
    MCoutput.Y = MCoutput.X;
end

end
