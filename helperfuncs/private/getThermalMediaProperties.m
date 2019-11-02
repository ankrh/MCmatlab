function model = getThermalMediaProperties(model)
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

mP_fH = model.HS.mediaProperties_funcHandles;

if ~isnan(model.HS.T(1))
  T = model.HS.T;
elseif ~isnan(model.HS.Tinitial(1))
  if numel(model.HS.Tinitial) == 1
    T = model.HS.Tinitial*ones(size(G.N_raw),'single');
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

%% Loop through different media and split if it has a dependence on T
M = zeros(size(M_raw),'uint8');
j = 1; % Position in mediaProperties
for i=1:length(mP_fH)
  paramVals = NaN(3,sum(M_raw(:) == i),'single'); % Number of columns is number of voxels with this medium, and rows are (VHC, YC, 1)
  if isa(mP_fH(i).VHC,'function_handle')
    paramVals(1,:) = mP_fH(i).VHC(T(M_raw(:) == i),FD(M_raw(:) == i));
  else
    paramVals(1,:) = mP_fH(i).VHC;
  end
  if isa(mP_fH(i).TC,'function_handle')
    paramVals(2,:) = mP_fH(i).TC(T(M_raw(:) == i),FD(M_raw(:) == i));
  else
    paramVals(2,:) = mP_fH(i).TC;
  end
  paramVals(3,:) = 1;
  
  if ~isa(mP_fH(i).VHC,'function_handle') && ...
     ~isa(mP_fH(i).TC,'function_handle') % No dependence
    M(M_raw(:) == i) = j;
    mP(j).name = mP_fH(i).name; % mP is mediaProperties
    mP(j).VHC = mP_fH(i).VHC;
    mP(j).TC = mP_fH(i).TC;
    mP(j).A = mP_fH(i).A;
    mP(j).E = mP_fH(i).E;
    j = j + 1;
  else
    if isnan(mP_fH(i).nBins)
      error('Error: nBins not specified for medium %s',mP_fH(i).name);
    end
    N = mP_fH(i).nBins; % Number of sub-media to split the i'th medium into
    [downsampledParamVals,binidx] = getDownsampledParamVals(paramVals,N);
    
    [~,sortidx] = sortrows(downsampledParamVals.');
    downsampledParamVals = downsampledParamVals(:,sortidx.'); % Sort the downsampled parameter values
    I_rev = zeros(N,1);
    I_rev(sortidx) = 1:N;
    binidx = I_rev(binidx); % binidx should be rearranged according to the same sort
    
    M(M_raw(:) == i) = j + binidx - 1;
    for k=N-1:-1:0 % Go backwards for slightly more efficient memory allocation (resizing only once instead of n_splitpts times)
      mP(j+k).name = mP_fH(i).name; % mP is mediaProperties
      mP(j+k).VHC = double(downsampledParamVals(1,k+1));
      mP(j+k).TC = double(downsampledParamVals(2,k+1));
      mP(j+k).A = mP_fH(i).A;
      mP(j+k).E = mP_fH(i).E;
    end
    j = j + N;
  end
end

if j-1 > 256
  error('Error: The total number of (sub-)media may not exceed 256');
end

%% Throw an error if a variable doesn't conform to its required interval
for j=1:length(mP)
  if mP(j).VHC <= 0
    error('Error: Medium %s has VHC <= 0',mP(j).name);
  elseif mP(j).TC < 0
    error('Error: Medium %s has TC < 0',mP(j).name);
  elseif mP(j).A < 0
    error('Error: Medium %s has A < 0',mP(j).name);
  elseif mP(j).E < 0
    error('Error: Medium %s has E < 0',mP(j).name);
  end
end

model.HS.mediaProperties = mP;
model.HS.M = M;
end