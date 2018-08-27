function runMonteCarloFluorescence(name)
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

%% Load data from defineGeometry.m and runMonteCarlo.m
load(['./Data/' name '.mat'],'G');
load(['./Data/' name '_MCoutput.mat'],'MCoutput');

%% USER SPECIFIED: Define simulation behavior
% Should the script run in silent mode? (disables overwrite prompt,
% command window text, progress indication and plot generation)
silentMode = false;

% Should MCmatlab use all available processors on Windows? Otherwise,
% one will be left unused. Useful for doing other work on the PC
% while simulations are running.
useAllCPUs = true;

% Simulation duration
simulationTime = 0.1;      % [min] time duration of the simulation

%% USER SPECIFIED: Define power of excitation beam
% Incident excitation power, for taking fluorescence saturations into account
P_excitation = 2; % [W]

%% USER SPECIFIED: Optional light collector properties
% A "light collector" in this context can be either an objective lens or a fiber tip
useLightCollector = false;

% Position of either the center of the objective lens focal plane or the fiber tip
xFPC_LC = 0; % [cm]
yFPC_LC = 0; % [cm]
zFPC_LC = 0.03; % [cm] negative values correspond to a location above the volume

% Direction that the light collector is facing, defined in the same way as the beam direction using
% ISO spherical coordinates.
% For above the surface (negative z), you'll want to satisfy 0<=theta<=pi/2.
% For below the surface (positive z), such as measuring in transmission, you'll want pi/2<=theta<=pi
% If theta = 0 or pi, phi serves only to rotate the view. If you want the light collector X axis 
% to coincide with the cuboid x axis, use phi = pi/2.
theta_LC = 0; % [rad]
phi_LC   = pi/2; % [rad]

% Focal length of the objective lens (if light collector is a fiber, set this to Inf).
f_LC = 1; % [cm]

% Diameter of the light collector aperture. For an ideal thin lens, this is 2*f*tan(asin(lensNA)).
diam_LC = 2; % [cm]

% Field Size of the imaging system (diameter of area in object plane that gets imaged). Only used for finite f_LC.
FieldSize_LC = 2; % [cm]

% Fiber NA. Only used for infinite f_LC.
NA_LC = 0.22; % [-]

% X and Y resolution of light collector in pixels
res_LC = 200;

%% Check for preexisting files
if(~silentMode && ~deleteDataFiles(name)); return; end

%% Calculate 3D fluorescence source distribution, including saturation
mua_vec = [G.mediaProperties.mua]; % The media's excitation absorption coefficients
Y_vec = [G.mediaProperties.Y]; % The media's fluorescence power yields
sat_vec = [G.mediaProperties.sat]; % The media's fluorescence saturation fluence rates (intensity)
sourceDistribution = Y_vec(G.M).*mua_vec(G.M)*P_excitation.*MCoutput.F./(1 + P_excitation*MCoutput.F./sat_vec(G.M)); % [W/cm^3]
if(max(sourceDistribution(:)) == 0); error('Error: No fluorescence emitters'); end

%% Check to ensure that the light collector is not inside the cuboid and set res_LC to 1 if using fiber
if useLightCollector
    if isfinite(f_LC)
        xLCC = xFPC_LC - f_LC*sin(theta_LC)*cos(phi_LC); % x position of Light Collector Center
        yLCC = yFPC_LC - f_LC*sin(theta_LC)*sin(phi_LC); % y position
        zLCC = zFPC_LC - f_LC*cos(theta_LC);             % z position
    else
        xLCC = xFPC_LC;
        yLCC = yFPC_LC;
        zLCC = zFPC_LC;
        res_LC = 1;
    end

    if (abs(xLCC)               < G.nx*G.dx/2 && ...
        abs(yLCC)               < G.ny*G.dy/2 && ...
        abs(zLCC - G.nz*G.dz/2) < G.nz*G.dz/2)
        error('Error: Light collector center (%.4f,%.4f,%.4f) is inside cuboid',xLCC,yLCC,zLCC);
    end    
end

%% Prepare structs
G.M = G.M - 1; % The medium matrix has to be converted from MATLAB's 1-based indexing to C's 0-based indexing
Beam = struct('P_excitation',P_excitation,'sourceDistribution',sourceDistribution);
LightCollector = struct('xFPC_LC',xFPC_LC,'yFPC_LC',yFPC_LC,'zFPC_LC',zFPC_LC,'theta_LC',theta_LC,'phi_LC',phi_LC,'f_LC',f_LC,...
    'diam_LC',diam_LC,'FieldSize_LC',FieldSize_LC,'NA_LC',NA_LC,'res_LC',res_LC);
MCinput_f = struct('silentMode',silentMode,'useAllCPUs',useAllCPUs,'simulationTime',simulationTime,...
    'useLightCollector',useLightCollector,'G',G,'Beam',Beam,'LightCollector',LightCollector);
clear G Beam LightCollector

%% Call Monte Carlo C script (MEX file) to get fluorescence fluence rate (intensity) distribution
MCoutput_f = MCmatlab(MCinput_f); % MCoutput_f.F is an absolute fluence rate (intensity) quantity, unlike the non-fluorescence MCoutput.F which are actually fluence rates normalized to the incident power
clear MCmatlab; % Unload MCmatlab MEX file so it can be modified externally again

%% Save output and clear memory
save(['./Data/' name '_MCoutput_fluorescence.mat'],'MCoutput_f','MCinput_f');
if(~silentMode); fprintf('./Data/%s_MCoutput_fluorescence.mat saved\n',name); end
clear MCinput_f MCoutput_f

%% Make plots
if(~silentMode); plotMCmatlab(name); end

end
