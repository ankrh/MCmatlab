function MCoutput = runMonteCarlo(mexMCinput)
%   Requires
%       MCmatlab.mex (architecture specific)
%
%   See also defineGeometry, plotMCmatlab, runMonteCarloFluorescence, simulateHeatDistribution

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

G = mexMCinput.G;

if ~isfield(mexMCinput,'wavelength')
  error('Error: No wavelength defined');
end
if ~isfield(mexMCinput,'matchedInterfaces')
  mexMCinput.matchedInterfaces = true;
end
if ~isfield(mexMCinput,'silentMode')
  mexMCinput.silentMode = false;
end
if ~isfield(mexMCinput,'useAllCPUs')
  mexMCinput.useAllCPUs = false;
end
if ~isfield(mexMCinput,'calcF')
  mexMCinput.calcF = true;
end
if ~isfield(mexMCinput,'calcFdet')
  mexMCinput.calcFdet = false;
end
if ~isfield(mexMCinput,'nExamplePaths')
  mexMCinput.nExamplePaths = 0;
end
if ~isfield(mexMCinput,'farfieldRes')
  mexMCinput.farfieldRes = 0;
end

mexMCinput.mediaProperties = getMediaProperties(G.mediaPropertiesFunc,mexMCinput.wavelength,G.mediaPropParams);

%% Extract the refractive indices if not assuming matched interfaces, otherwise assume all 1's
if(~mexMCinput.matchedInterfaces)
  for j=1:length(mexMCinput.mediaProperties) % Check that all media have a refractive index defined
    if(~isfield(mexMCinput.mediaProperties,'n') || any(isempty(mexMCinput.mediaProperties(j).n)))
      error('matchedInterfaces is false, but refractive index isn''t defined for all media');
    end
  end
  n_vec = [mexMCinput.mediaProperties.n];
  for j=1:G.nz % Check that each xy slice has constant refractive index, so refractive index is only a function of z
    if(length(unique(n_vec(G.M(:,:,j)))) > 1)
      error('matchedInterfaces is false, but refractive index isn''t constant for z index %d (z = %f).\nEach xy slice must have constant refractive index.',j,G.z(j));
    end
  end
  mexMCinput.RI = n_vec(G.M(1,1,:));
else
  [mexMCinput.mediaProperties.n] = deal(1);
  mexMCinput.RI = ones(G.nz,1);
end

if isfield(mexMCinput,'LightCollector')
  %% Check to ensure that the light collector is not inside the cuboid and set res to 1 if using fiber
  mexMCinput.useLightCollector = true;
  if mexMCinput.boundaryType == 0
    error('Error: If boundaryType == 0, no photons can escape to be registered on the light collector. Disable light collector or change boundaryType.');
  end
  if isfinite(mexMCinput.LightCollector.f)
    xLCC = mexMCinput.LightCollector.x - mexMCinput.LightCollector.f*sin(mexMCinput.LightCollector.theta)*cos(mexMCinput.LightCollector.phi); % x position of Light Collector Center
    yLCC = mexMCinput.LightCollector.y - mexMCinput.LightCollector.f*sin(mexMCinput.LightCollector.theta)*sin(mexMCinput.LightCollector.phi); % y position
    zLCC = mexMCinput.LightCollector.z - mexMCinput.LightCollector.f*cos(mexMCinput.LightCollector.theta);                                 % z position
  else
    xLCC = mexMCinput.LightCollector.x;
    yLCC = mexMCinput.LightCollector.y;
    zLCC = mexMCinput.LightCollector.z;
    if mexMCinput.LightCollector.res ~= 1
      error('Error: LightCollector.res must be 1 when LightCollector.f is Inf');
    end
  end

  if (abs(xLCC)               < G.nx*G.dx/2 && ...
      abs(yLCC)               < G.ny*G.dy/2 && ...
      abs(zLCC - G.nz*G.dz/2) < G.nz*G.dz/2)
    error('Error: Light collector center (%.4f,%.4f,%.4f) is inside cuboid',xLCC,yLCC,zLCC);
  end

  %% If no time tagging bins are defined, assume no time tagging is to be performed
  if ~isfield(mexMCinput.LightCollector,'nTimeBins')
    mexMCinput.LightCollector.tStart = 0;
    mexMCinput.LightCollector.tEnd = 0;
    mexMCinput.LightCollector.nTimeBins = 0;
  end
else
  %% Assume no light collector is present
  mexMCinput.useLightCollector = false;
  mexMCinput.LightCollector.x = 0;
  mexMCinput.LightCollector.y = 0;
  mexMCinput.LightCollector.z = 0;
  mexMCinput.LightCollector.theta = 0;
  mexMCinput.LightCollector.phi = 0;
  mexMCinput.LightCollector.f = 0;
  mexMCinput.LightCollector.diam = 0;
  mexMCinput.LightCollector.FieldSize = 0;
  mexMCinput.LightCollector.NA = 0;
  mexMCinput.LightCollector.res = 0;
  mexMCinput.LightCollector.tStart = 0;
  mexMCinput.LightCollector.tEnd = 0;
  mexMCinput.LightCollector.nTimeBins = 0;
end

if xor(isfield(mexMCinput.Beam,'nearFieldType'), isfield(mexMCinput.Beam,'farFieldType'))
  error('Error: nearFieldType and farFieldType must either both be specified, or neither');
end
if isfield(mexMCinput,'simulationTime') && isfield(mexMCinput,'nPhotons')
  error('Error: simulationTime and nPhotons may not both be specified');
end
if isfield(mexMCinput.Beam,'nearFieldType') && isfield(mexMCinput.Beam,'beamType')
  error('Error: nearFieldType and beamType may not both be specified');
end
if ~mexMCinput.calcF && ~mexMCinput.useLightCollector
  error('Error: calcF is false, but no light collector is defined');
end
if mexMCinput.calcFdet && ~mexMCinput.useLightCollector
  error('Error: calcFdet is true, but no light collector is defined');
end
if mexMCinput.farfieldRes && mexMCinput.boundaryType == 0
  error('Error: If boundaryType == 0, no photons can escape to be registered in the far field. Set farfieldRes to zero or change boundaryType.');
end

if isfield(mexMCinput.Beam,'nearFieldType')
  mexMCinput.Beam.beamType = -1;
else
  mexMCinput.Beam.nearFieldType = -1;
  mexMCinput.Beam.farFieldType = -1;
end

%% Call Monte Carlo C script (MEX file) to get fluence rate (intensity) distribution
mexMCinput.G.M = mexMCinput.G.M-1; % Convert to C-style indexing
MCoutput = MCmatlab(mexMCinput);

% Add positions of the centers of the pixels in the light collector image
if mexMCinput.useLightCollector && mexMCinput.LightCollector.res > 1
  MCoutput.X = linspace(mexMCinput.LightCollector.FieldSize*(1/mexMCinput.LightCollector.res-1),mexMCinput.LightCollector.FieldSize*(1-1/mexMCinput.LightCollector.res),mexMCinput.LightCollector.res)/2;
  MCoutput.Y = MCoutput.X;
end

% Add angles of the centers of the far field pixels
if mexMCinput.farfieldRes
  MCoutput.FarFieldTheta = linspace(pi/mexMCinput.farfieldRes/2,pi-pi/mexMCinput.farfieldRes/2,mexMCinput.farfieldRes);
  MCoutput.FarFieldPhi   = linspace(-pi+(2*pi)/mexMCinput.farfieldRes/2,pi-(2*pi)/mexMCinput.farfieldRes/2,mexMCinput.farfieldRes);
end

% Also put the mediaProperties into the MCoutput struct for later use by other functions
MCoutput.mediaProperties = mexMCinput.mediaProperties;
end
