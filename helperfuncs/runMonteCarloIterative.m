function MCoutput = runMonteCarloIterative(mexMCinput)
%   Requires
%       runMonteCarlo.m
%
%   See also defineGeometry, plotMCmatlab, runMonteCarlo, runMonteCarloFluorescence, simulateHeatDistribution

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

mexMCinput.mediaProperties_funcHandles = getmediaProperties_funcHandles(mexMCinput);

mexMCinput.mediaProperties = getmediaProperties(mexMCinput);

% Fill in fluorescence and Arrhenius parameter assumptions
% For all media for which the fluorescence power yield Y, Arrhenius
% activation energy E or Arrhenius pre-exponential factor A was not 
% specified, assume they are zero.
for j=1:length(mediaProperties_3)
  if(~isfield(mediaProperties_3,'Y') || isempty(mediaProperties_3(j).Y))
    mediaProperties_3(j).Y = 0;
  end
  if(~isfield(mediaProperties_3,'E') || isempty(mediaProperties_3(j).E))
    mediaProperties_3(j).E = 0;
  end
  if(~isfield(mediaProperties_3,'A') || isempty(mediaProperties_3(j).A))
    mediaProperties_3(j).A = 0;
  end
end

% Throw an error if a variable doesn't conform to its required interval
for j=1:length(mediaProperties_3)
  if(~isfield(mediaProperties_3,'mua') || isempty(mediaProperties_3(j).mua))
    error('Medium %s has no mua.',mediaProperties_3(j).name);
  elseif(~isfield(mediaProperties_3,'mus') || isempty(mediaProperties_3(j).mus))
    error('Medium %s has no mus.',mediaProperties_3(j).name);
  elseif(~isfield(mediaProperties_3,'g') || isempty(mediaProperties_3(j).g))
    error('Medium %s has no g.',mediaProperties_3(j).name);
  end

  if(mediaProperties_3(j).mua <= 0)
    error('Medium %s has mua <= 0',mediaProperties_3(j).name);
  elseif(mediaProperties_3(j).mus <= 0)
    error('Medium %s has mus <= 0',mediaProperties_3(j).name);
  elseif(abs(mediaProperties_3(j).g) > 1)
    error('Medium %s has abs(g) > 1',mediaProperties_3(j).name);
  elseif(isfield(mediaProperties_3,'n') && ~isempty(mediaProperties_3(j).n) && mediaProperties_3(j).n < 1)
    error('Medium %s has n < 1',mediaProperties_3(j).name);
  elseif(isfield(mediaProperties_3,'VHC') && ~isempty(mediaProperties_3(j).VHC) && mediaProperties_3(j).VHC <= 0)
    error('Medium %s has VHC <= 0',mediaProperties_3(j).name);
  elseif(isfield(mediaProperties_3,'TC') && ~isempty(mediaProperties_3(j).TC) && mediaProperties_3(j).TC < 0)
    error('Medium %s has TC < 0',mediaProperties_3(j).name);
  elseif(mediaProperties_3(j).Y < 0)
    error('Medium %s has Y < 0',mediaProperties_3(j).name);
  elseif(mediaProperties_3(j).E < 0)
    error('Medium %s has E < 0',mediaProperties_3(j).name);
  elseif(mediaProperties_3(j).A < 0)
    error('Medium %s has A < 0',mediaProperties_3(j).name);
  end
end

% Extract the refractive indices if not assuming matched interfaces, otherwise assume all 1's
if(~mexMCinput.matchedInterfaces)
  for j=1:length(mexMCinput.mediaProperties_3) % Check that all media have a refractive index defined
    if(~isfield(mexMCinput.mediaProperties_3,'n') || any(isempty(mexMCinput.mediaProperties_3(j).n)))
      error('matchedInterfaces is false, but refractive index isn''t defined for all media');
    end
  end
  n_vec = [mexMCinput.mediaProperties_3.n];
  for j=1:G.nz % Check that each xy slice has constant refractive index, so refractive index is only a function of z
    if(length(unique(n_vec(G.M(:,:,j)))) > 1)
      error('matchedInterfaces is false, but refractive index isn''t constant for z index %d (z = %f).\nEach xy slice must have constant refractive index.',j,G.z(j));
    end
  end
  mexMCinput.RI = n_vec(G.M(1,1,:));
else
  [mexMCinput.mediaProperties_3.n] = deal(1);
  mexMCinput.RI = ones(G.nz,1);
end

% Put media properties into the mexMCinput struct
mexMCinput.mediaProperties = mediaProperties_3;

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
