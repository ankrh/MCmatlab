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

mexMCinput = checkMCinputFields(mexMCinput);
mexMCinput = getMediaProperties_funcHandles(mexMCinput); % Calls the mediaPropertiesFunc and converts those fields that are specified as cell array formulas into function handles taking intensity and temperature as input
if mexMCinput.temperatureDependent && ~isfield(mexMCinput,'Tinitial')
  error('Error: Parameters are temperature dependent but no initial temperature is provided');
else
  mexMCinput.Tinitial = zeros(size(G.M)); % If there's no temperature dependence then it doesn't matter what we set these values to
end
if ~isfield(mexMCinput,'Iinitial')
  mexMCinput.Iinitial = zeros(size(G.M));
end
[M_split, mexMCinput.mediaProperties] = getMediaProperties(G.M,mexMCinput.mediaProperties_funcHandles,mexMCinput.Iinitial,mexMCinput.Tinitial);


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
