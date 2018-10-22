function FMCoutput = runMonteCarloFluorescence(FMCinput)
%   Script for simulating distribution and magnitude of fluorescence
%   based on the output of runMonteCarlo.m
%
%   Prepares and runs the Monte Carlo simulation.
%
%   Requires
%       MCmatlab.mex (architecture specific)
%
%	See also runMonteCarlo, plotMCmatlabFluorescence

%%%%%
%   Copyright 2018 by Anders K. Hansen, DTU Fotonik
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

if ~isfield(FMCinput,'silentMode')
	FMCinput.silentMode = false;
end
if ~isfield(FMCinput,'useAllCPUs')
	FMCinput.useAllCPUs = false;
end

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

%% Call Monte Carlo C script (MEX file) to get fluence rate (intensity) distribution
FMCinput.G.M = FMCinput.G.M - 1; % Convert to C-style indexing
FMCoutput = MCmatlab(FMCinput); % FMCoutput.F is an absolute fluence rate (intensity) quantity, unlike the non-fluorescence MCoutput.F which are actually fluence rates normalized to the incident power

% Add positions of the centers of the pixels in the light collector image
if FMCinput.useLightCollector && FMCinput.LightCollector.res > 1
	FMCoutput.X = linspace(FMCinput.LightCollector.FieldSize*(1/FMCinput.LightCollector.res-1),FMCinput.LightCollector.FieldSize*(1-1/FMCinput.LightCollector.res),FMCinput.LightCollector.res)/2;
    FMCoutput.Y = FMCoutput.X;
end

end
