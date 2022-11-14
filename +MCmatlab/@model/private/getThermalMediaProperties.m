function model = getThermalMediaProperties(model)
  %   Returns the reduced medium matrix, using only numbers from 1 up to the number of used media, and
  %   the known media properties (optical, thermal and/or fluorescence) at the specified wavelength.
  %
  %   See also mediaPropertiesLibrary

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
  mP_fH = G.mediaPropertiesFunc(G.mediaPropParams);
  [uniqueMedia,~,M_trim] = unique(G.M_raw);
  M_trim = reshape(M_trim,model.G.nx,model.G.ny,model.G.nz);
  mP_fHtrim = mP_fH(uniqueMedia);
  nM = numel(mP_fHtrim);
  if numel(unique({mP_fHtrim.name})) ~= nM
    error('Error: Two media have the same name.');
  end

  T = model.HS.T;

  FD = NaN([G.nx G.ny G.nz],'single');
  FD(:) = 1 - exp(-model.HS.Omega); % Fractional damage of molecules/cells. model.HS.Omega may be scalar 0 or 3D

  %% Loop through different media and split if it has a dependence on T
  nSMtotal = 0;
  for iM = 1:nM
    nSMtotal = nSMtotal + mP_fHtrim(iM).nBins;
  end
  if nSMtotal > 256
    error('Error: The total number of (sub-)media may not exceed 256');
  end
  mP = struct();
  mP.name = cell(nSMtotal,1);
  mP.VHC  = NaN(nSMtotal,1);
  mP.TC   = NaN(nSMtotal,1);
  mP.A    = NaN(nSMtotal,1);
  mP.E    = NaN(nSMtotal,1);

  VHC3D = NaN(G.nx, G.ny, G.nz, 'single');
  TC3D = NaN(G.nx, G.ny, G.nz, 'single');
  M = zeros(size(M_trim),'uint8');
  iSM = 1; % Position in mediaProperties
  for iM=1:nM
    idxs3D = M_trim(:) == iM;
    try VHC3D(idxs3D) = mP_fHtrim(iM).VHC(T(idxs3D),FD(idxs3D));
    catch; error('Error: The VHC function of medium %s is throwing an error. If you have specified the function with two input arguments, it must be able to accept both arguments (temperature T and fractional damage FD) as 3D arrays. Therefore, remember to use element-wise arithmetic operators (.* instead of *, ./ instead of / etc.)',mP_fHtrim(iM).name);
    end
    try TC3D(idxs3D) = mP_fHtrim(iM).TC(T(idxs3D),FD(idxs3D));
    catch; error('Error: The TC function of medium %s is throwing an error. If you have specified the function with two input arguments, it must be able to accept both arguments (temperature T and fractional damage FD) as 3D arrays. Therefore, remember to use element-wise arithmetic operators (.* instead of *, ./ instead of / etc.)',mP_fHtrim(iM).name);
    end
    if any(~isfinite(VHC3D(idxs3D))) || any(VHC3D(:) <= 0)
      error('The VHC function of medium %s returns values that are not finite and positive.',mP_fHtrim(iM).name);
    elseif any(~isfinite(TC3D(idxs3D))) || any(TC3D(:) <= 0)
      error('The TC function of medium %s returns values that are not finite and positive.',mP_fHtrim(iM).name);
    end

    paramVals = NaN(3,sum(M_trim(:) == iM),'single'); % Number of columns is number of voxels with this medium, and rows are (VHC, YC, 1)
    paramVals(1,:) = VHC3D(idxs3D);
    paramVals(2,:) = TC3D(idxs3D);
    paramVals(3,:) = 1;
    nSM = mP_fHtrim(iM).nBins; % Number of sub-media to split the iM'th medium into
    [downsampledParamVals,binidx] = getDownsampledParamVals(paramVals,nSM);
    [~,sortidx] = sortrows(downsampledParamVals.');
    downsampledParamVals = downsampledParamVals(:,sortidx.'); % Sort the downsampled parameter values
    I_rev = zeros(nSM,1);
    I_rev(sortidx) = 1:nSM;
    binidx = I_rev(binidx); % binidx is rearranged according to the same sort

    M(idxs3D) = iSM + binidx - 1;
    mP.name(iSM:iSM+nSM-1) = {deal(mP_fHtrim(iM).name)};
    mP.VHC (iSM:iSM+nSM-1) = double(downsampledParamVals(1,:));
    mP.TC  (iSM:iSM+nSM-1) = double(downsampledParamVals(2,:));
    mP.A   (iSM:iSM+nSM-1) = deal(mP_fHtrim(iM).A);
    mP.E   (iSM:iSM+nSM-1) = deal(mP_fHtrim(iM).E);
    iSM = iSM + nSM;
  end

  model.HS.mediaProperties = mP;
  model.HS.M = M;
end