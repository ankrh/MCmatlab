function FMCoutput = runMonteCarloFluorescence(mexFMCinput)
%   Script for simulating distribution and magnitude of fluorescence
%   based on the output of runMonteCarlo.m
%
%   Prepares and runs the Monte Carlo simulation.
%
%   Requires
%       MCmatlab.mex (architecture specific)
%
%   See also runMonteCarlo, plotMCmatlabFluorescence

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

G = mexFMCinput.G;

if ~isfield(mexFMCinput,'wavelength')
  error('Error: No wavelength defined');
end
if ~isfield(mexFMCinput,'silentMode')
  mexFMCinput.silentMode = false;
end
if ~isfield(mexFMCinput,'useAllCPUs')
  mexFMCinput.useAllCPUs = false;
end
if ~isfield(mexFMCinput,'calcF')
  mexFMCinput.calcF = true;
end
if ~isfield(mexFMCinput.MCoutput,'F')
  error('Error: F matrix not calculated for excitation light');
end
if ~isfield(mexFMCinput,'calcFdet')
  mexFMCinput.calcFdet = false;
end
if ~isfield(mexFMCinput,'nExamplePaths')
  mexFMCinput.nExamplePaths = 0;
end
if ~isfield(mexFMCinput,'farfieldRes')
  mexFMCinput.farfieldRes = 0;
end

mexFMCinput.mediaProperties = getMediaProperties(G.mediaPropertiesFunc,mexFMCinput.wavelength,G.mediaPropParams);

%% Extract the refractive indices if not assuming matched interfaces, otherwise assume all 1's
if(~mexFMCinput.matchedInterfaces)
  for j=1:length(mexFMCinput.mediaProperties) % Check that all media have a refractive index defined
    if(~isfield(mexFMCinput.mediaProperties,'n') || any(isempty(mexFMCinput.mediaProperties(j).n)))
      error('matchedInterfaces is false, but refractive index isn''t defined for all media');
    end
  end
  n_vec = [mexFMCinput.mediaProperties.n];
  for j=1:G.nz % Check that each xy slice has constant refractive index, so refractive index is only a function of z
    if(length(unique(n_vec(G.M(:,:,j)))) > 1)
      error('matchedInterfaces is false, but refractive index isn''t constant for z index %d (z = %f).\nEach xy slice must have constant refractive index.',j,G.z(j));
    end
  end
  mexFMCinput.RI = n_vec(G.M(1,1,:));
else
  [mexFMCinput.mediaProperties.n] = deal(1);
  mexFMCinput.RI = ones(G.nz,1);
end

%% Calculate 3D fluorescence source distribution
mua_vec = [mexFMCinput.MCoutput.mediaProperties.mua]; % The media's excitation absorption coefficients
Y_vec = [mexFMCinput.MCoutput.mediaProperties.Y]; % The media's fluorescence power yields
mexFMCinput.Beam.sourceDistribution = Y_vec(G.M).*mua_vec(G.M).*mexFMCinput.MCoutput.F; % [W/cm^3]
if(max(mexFMCinput.Beam.sourceDistribution(:)) == 0); error('Error: No fluorescence emitters'); end

%% Check to ensure that the light collector is not inside the cuboid and set res_LC to 1 if using fiber
if isfield(mexFMCinput,'LightCollector')
  mexFMCinput.useLightCollector = true;
  if mexFMCinput.boundaryType == 0
    error('Error: If boundaryType == 0, no photons can escape to be registered on the light collector. Disable light collector or change boundaryType.');
  end
  if isfinite(mexFMCinput.LightCollector.f)
    xLCC = mexFMCinput.LightCollector.x - mexFMCinput.LightCollector.f*sin(mexFMCinput.LightCollector.theta)*cos(mexFMCinput.LightCollector.phi); % x position of Light Collector Center
    yLCC = mexFMCinput.LightCollector.y - mexFMCinput.LightCollector.f*sin(mexFMCinput.LightCollector.theta)*sin(mexFMCinput.LightCollector.phi); % y position
    zLCC = mexFMCinput.LightCollector.z - mexFMCinput.LightCollector.f*cos(mexFMCinput.LightCollector.theta);                                  % z position
  else
    xLCC = mexFMCinput.LightCollector.x;
    yLCC = mexFMCinput.LightCollector.y;
    zLCC = mexFMCinput.LightCollector.z;
    if mexFMCinput.LightCollector.res ~= 1
      error('Error: LightCollector.res must be 1 when LightCollector.f is Inf');
    end
  end

  if (abs(xLCC)               < G.nx*G.dx/2 && ...
      abs(yLCC)               < G.ny*G.dy/2 && ...
      abs(zLCC - G.nz*G.dz/2) < G.nz*G.dz/2)
    error('Error: Light collector center (%.4f,%.4f,%.4f) is inside cuboid',xLCC,yLCC,zLCC);
  end

  if ~isfield(mexFMCinput.LightCollector,'nTimeBins')
    mexFMCinput.LightCollector.tStart = 0;
    mexFMCinput.LightCollector.tEnd = 0;
    mexFMCinput.LightCollector.nTimeBins = 0;
  elseif mexFMCinput.LightCollector.nTimeBins > 0
    error('Error: Fluorescence Monte Carlo does not support time tagging: LightCollector.nTimeBins must be 0.');
  end
else
  mexFMCinput.useLightCollector = false;
  mexFMCinput.LightCollector.x = 0;
  mexFMCinput.LightCollector.y = 0;
  mexFMCinput.LightCollector.z = 0;
  mexFMCinput.LightCollector.theta = 0;
  mexFMCinput.LightCollector.phi = 0;
  mexFMCinput.LightCollector.f = 0;
  mexFMCinput.LightCollector.diam = 0;
  mexFMCinput.LightCollector.FieldSize = 0;
  mexFMCinput.LightCollector.NA = 0;
  mexFMCinput.LightCollector.res = 0;
  mexFMCinput.LightCollector.tStart = 0;
  mexFMCinput.LightCollector.tEnd = 0;
  mexFMCinput.LightCollector.nTimeBins = 0;
end

if mexFMCinput.farfieldRes && mexFMCinput.boundaryType == 0
  error('Error: If boundaryType == 0, no photons can escape to be registered in the far field. Set farfieldRes to zero or change boundaryType.');
end
if ~mexFMCinput.calcF && ~mexFMCinput.useLightCollector
  error('Error: calcF is false, but no light collector is defined');
end
if mexFMCinput.calcFdet && ~mexFMCinput.useLightCollector
  error('Error: calcFdet is true, but no light collector is defined');
end

%% Call Monte Carlo C script (MEX file) to get fluence rate (intensity) distribution
mexFMCinput.G.M = mexFMCinput.G.M - 1; % Convert to C-style indexing
FMCoutput = MCmatlab(mexFMCinput); % FMCoutput.F is the fluence rate normalized to the incident excitation power (not the emitted fluorescence power)

% Add positions of the centers of the pixels in the light collector image
if mexFMCinput.useLightCollector && mexFMCinput.LightCollector.res > 1
  FMCoutput.X = linspace(mexFMCinput.LightCollector.FieldSize*(1/mexFMCinput.LightCollector.res-1),mexFMCinput.LightCollector.FieldSize*(1-1/mexFMCinput.LightCollector.res),mexFMCinput.LightCollector.res)/2;
  FMCoutput.Y = FMCoutput.X;
end

% Add angles of the centers of the far field pixels
if mexFMCinput.farfieldRes
  FMCoutput.FarFieldTheta = linspace(pi/mexFMCinput.farfieldRes/2,pi-pi/mexFMCinput.farfieldRes/2,mexFMCinput.farfieldRes);
  FMCoutput.FarFieldPhi   = linspace(-pi+(2*pi)/mexFMCinput.farfieldRes/2,pi-(2*pi)/mexFMCinput.farfieldRes/2,mexFMCinput.farfieldRes);
end

% Also put the mediaProperties into the FMCoutput struct for later use by other functions
FMCoutput.mediaProperties = mexFMCinput.mediaProperties;
end
