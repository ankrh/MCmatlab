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

%% Check to ensure that the light collector is not inside the cuboid and set res_LC to 1 if using fiber
if MCinput.useLightCollector
    if isfinite(MCinput.LightCollector.f_LC)
        xLCC = MCinput.LightCollector.xFPC_LC - MCinput.LightCollector.f_LC*sin(MCinput.LightCollector.theta_LC)*cos(MCinput.LightCollector.phi_LC); % x position of Light Collector Center
        yLCC = MCinput.LightCollector.yFPC_LC - MCinput.LightCollector.f_LC*sin(MCinput.LightCollector.theta_LC)*sin(MCinput.LightCollector.phi_LC); % y position
        zLCC = MCinput.LightCollector.zFPC_LC - MCinput.LightCollector.f_LC*cos(MCinput.LightCollector.theta_LC);             % z position
    else
        xLCC = MCinput.LightCollector.xFPC_LC;
        yLCC = MCinput.LightCollector.yFPC_LC;
        zLCC = MCinput.LightCollector.zFPC_LC;
        MCinput.LightCollector.res_LC = 1;
    end

    if (abs(xLCC)                           < MCinput.G.nx*MCinput.G.dx/2 && ...
        abs(yLCC)                           < MCinput.G.ny*MCinput.G.dy/2 && ...
        abs(zLCC - MCinput.G.nz*MCinput.G.dz/2) < MCinput.G.nz*MCinput.G.dz/2)
        error('Error: Light collector center (%.4f,%.4f,%.4f) is inside cuboid',xLCC,yLCC,zLCC);
    end    
end

%% Call Monte Carlo C script (MEX file) to get fluence rate (intensity) distribution
MCinput.G.M = MCinput.G.M-1; % Convert to C-style indexing
MCoutput = MCmatlab(MCinput);
clear MCmatlab; % Unload MCmatlab MEX file so it can be modified externally again
end
