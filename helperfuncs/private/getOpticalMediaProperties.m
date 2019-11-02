function model = getOpticalMediaProperties(model,simType)
%   Returns the reduced medium matrix, using only numbers from 1 up to the number of used media, and
%   the known media properties (optical, thermal and/or fluorescence) at the specified wavelength.
%
%   See also mediaPropertiesLibrary, defineGeometry

%%%%%
%   Copyright 2017, 2018 by Dominik Marti and Anders K. Hansen, DTU Fotonik
%   This function was inspired by makeTissueList.m of the mcxyz program hosted at omlc.org
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

G = model.G;

M_raw = G.M_raw;

if simType == 2
  mP_fH = model.FMC.mediaProperties_funcHandles;
  matchedInterfaces = model.FMC.matchedInterfaces;
else
  mP_fH = model.MC.mediaProperties_funcHandles;
  matchedInterfaces = model.MC.matchedInterfaces;
end

if isnan(model.MC.FR(1))
  FR = zeros(size(G.M_raw),'single');
else
  FR = model.MC.FR;
end
if ~isnan(model.HS.T(1))
  T = model.HS.T;
elseif ~isnan(model.HS.Tinitial(1))
  if numel(model.HS.Tinitial) == 1
    T = model.HS.Tinitial*ones(size(G.M_raw),'single');
  else
    T = single(model.HS.Tinitial);
  end
else
  T = zeros(size(G.M_raw),'single');
end

if ~isnan(model.HS.Omega(1))
  FD = 1 - exp(-model.HS.Omega); % Fractional damage of molecules/cells
else
  FD = zeros(size(G.M_raw),'single');
end

%% Loop through different media and split if it has a dependence on FR or T
M = zeros(size(M_raw),'uint8');
j = 1; % Position in mediaProperties
for i=1:length(mP_fH)
  paramVals = NaN(3,sum(M_raw(:) == i),'single'); % Number of columns is number of voxels with this medium, and rows are (mua, mus, g)
  if isa(mP_fH(i).mua,'function_handle')
    paramVals(1,:) = mP_fH(i).mua(FR(M_raw(:) == i),T(M_raw(:) == i),FD(M_raw(:) == i));
  else
    paramVals(1,:) = mP_fH(i).mua;
  end
  if isa(mP_fH(i).mus,'function_handle')
    paramVals(2,:) = mP_fH(i).mus(FR(M_raw(:) == i),T(M_raw(:) == i),FD(M_raw(:) == i));
  else
    paramVals(2,:) = mP_fH(i).mus;
  end
  if isa(mP_fH(i).g,'function_handle')
    paramVals(3,:) = mP_fH(i).g(FR(M_raw(:) == i),T(M_raw(:) == i),FD(M_raw(:) == i));
  else
    paramVals(3,:) = mP_fH(i).g;
  end
  
  if ~isa(mP_fH(i).mua,'function_handle') && ...
     ~isa(mP_fH(i).mus,'function_handle') && ...
     ~isa(mP_fH(i).g  ,'function_handle') % No dependence
    M(M_raw(:) == i) = j;
    mP(j).name = mP_fH(i).name; % mP is mediaProperties
    mP(j).mua = mP_fH(i).mua;
    mP(j).mus = mP_fH(i).mus;
    mP(j).g = mP_fH(i).g;
    mP(j).n = mP_fH(i).n;
    mP(j).Y = mP_fH(i).Y; % May still be a function handle
    j = j + 1;
  else
    if isnan(mP_fH(i).nBins)
      error('Error: nBins not specified for medium %s',mP_fH(i).name);
    end
    N = mP_fH(i).nBins; % Number of sub-media to split the i'th medium into
    [downsampledParamVals,binidx] = getDownsampledParamVals(paramVals,N); % paramVals and downsampledParamVals are single precision and binidx is uint8
    
    [~,sortidx] = sortrows(downsampledParamVals.');
    downsampledParamVals = downsampledParamVals(:,sortidx.'); % Sort the downsampled parameter values
    I_rev = zeros(N,1);
    I_rev(sortidx) = 1:N;
    binidx = I_rev(binidx); % binidx should be rearranged according to the same sort
    
    M(M_raw(:) == i) = j + binidx - 1;
    for k=N-1:-1:0 % Go backwards for slightly more efficient memory allocation (resizing only once instead of n_splitpts times)
      mP(j+k).name = mP_fH(i).name; % mP is mediaProperties
      mP(j+k).mua = double(downsampledParamVals(1,k+1));
      mP(j+k).mus = double(downsampledParamVals(2,k+1));
      mP(j+k).g   = double(downsampledParamVals(3,k+1));
      mP(j+k).n = mP_fH(i).n;
      mP(j+k).Y = mP_fH(i).Y; % May still be a function handle
    end
    j = j + N;
  end
end

if j-1 > 256
  error('Error: The total number of (sub-)media may not exceed 256');
end

%% Throw an error if a variable doesn't conform to its required interval
for j=1:length(mP)
  if ~isfinite(mP(j).mua) || mP(j).mua <= 0
    error('Error: Medium %s has invalid mua (%f)',mP(j).name,mP(j).mua);
  elseif ~isfinite(mP(j).mus) || mP(j).mus <= 0
    error('Error: Medium %s has invalid mus (%f)',mP(j).name,mP(j).mus);
  elseif ~isfinite(mP(j).g) || abs(mP(j).g) > 1
    error('Error: Medium %s has invalid g (%f)',mP(j).name,mP(j).g);
  elseif mP(j).n < 1
    error('Error: Medium %s has invalid n (%f)',mP(j).name,mP(j).n);
  end
end

%% Extract the refractive indices if not assuming matched interfaces, otherwise assume all 1's
if ~matchedInterfaces
  n_vec = [mP.n];
  if any(isnan(n_vec)) % Check that all media have a refractive index defined
    error('Error: matchedInterfaces is false, but refractive index isn''t defined for all media');
  end
  for j=1:G.nz % Check that each xy slice has constant refractive index, so refractive index is only a function of z
    if length(unique(n_vec(M(:,:,j)))) > 1
      error('Error: matchedInterfaces is false, but refractive index isn''t constant for z index %d (z = %f).\nEach xy slice must have constant refractive index.',j,G.z(j));
    end
  end
  RI = n_vec(M(1,1,:));
else
  RI = ones(G.nz,1);
end

if simType == 1
  model.MC.mediaProperties = mP;
  model.MC.RI = RI;
  model.MC.M = M;
else
  model.FMC.mediaProperties = mP;
  model.FMC.RI = RI;
  model.FMC.M = M;
end
end