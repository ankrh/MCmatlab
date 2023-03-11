function model = getOpticalMediaProperties(model,simType)
  G = model.G;
  mP_fH = G.mediaPropertiesFunc(G.mediaPropParams);
  [uniqueMedia,~,M_trim] = unique(G.M_raw);
  M_trim = reshape(M_trim,model.G.nx,model.G.ny,model.G.nz);
  mP_fHtrim = mP_fH(uniqueMedia);
  nM = numel(mP_fHtrim);
  if numel(unique({mP_fHtrim.name})) ~= nM
    error('Error: Two media have the same name.');
  end

  if simType == 2
    matchedInterfaces = model.FMC.matchedInterfaces;
    smoothingLengthScale = model.FMC.smoothingLengthScale;
    Ls = model.FMC.wavelength; % lambdas
    DC = model.FMC.depositionCriteria;
  else
    matchedInterfaces = model.MC.matchedInterfaces;
    smoothingLengthScale = model.MC.smoothingLengthScale;
    Ls = model.MC.wavelength;
    DC = model.MC.depositionCriteria;
  end

  T = NaN([G.nx G.ny G.nz],'single');
  T(:) = model.HS.T; % model.HS.T may be scalar or 3D

  FD = NaN([G.nx G.ny G.nz],'single');
  FD(:) = 1 - exp(-model.HS.Omega); % Fractional damage of molecules/cells. model.HS.Omega may be scalar 0 or 3D

  FR = NaN([G.nx G.ny G.nz],'single');
  FR(:) = model.MC.P*sum(model.MC.NFR,4); % All excitation wavelengths are summed to get the total fluence rate. model.MC.NFR may be scalar 0 or 3D or 4D

  %% Loop through wavelengths
  uniqueCDFs = [];
  M = zeros(size(M_trim),'uint8');
  nL = numel(Ls);
  nSMtotal = 0;
  for iM = 1:nM
    nSMtotal = nSMtotal + mP_fHtrim(iM).nBins;
  end
  if nSMtotal > 256
    error('Error: The total number of (sub-)media may not exceed 256');
  end
  mP = struct();
  mP.name   = cell(nSMtotal,1);
  mP.mua    = NaN(nSMtotal,nL);
  mP.mus    = NaN(nSMtotal,nL);
  mP.g      = NaN(nSMtotal,nL);
  mP.n      = NaN(nSMtotal,nL);
  mP.CDFidx = NaN(nSMtotal,nL);

  %% Loop through different media and split if it has a dependence on FR or T
  iSM = 1; % Index of sub-medium, the first dimension position in the output mediaProperties
  for iM=1:nM % Index of medium, the position in the trimmed input mediaProperties
    idxs3D = M_trim(:) == iM;
    if any(G.FRdependent(iM,1:3)) || any(G.optTdependent(iM,1:3)) || any(G.optFDdependent(iM,1:3))
      FRs = FR(idxs3D);
      Ts = T(idxs3D);
      FDs = FD(idxs3D);
    end
    for iL = 1:nL % The second dimension position in the output mediaProperties
      L = Ls(iL);
      try
        if G.FRdependent(iM,1) || G.optTdependent(iM,1) || G.optFDdependent(iM,1)
          muavals = single(mP_fHtrim(iM).mua(L,FRs,Ts,FDs));
        else
          muavals = single(mP_fHtrim(iM).mua(L,0,0,0));
        end
      catch; error('Error: The mua function of medium %s is throwing an error. If you have specified the function with four input arguments, it must be able to accept the first (wavelength in nm) as a scalar and the next three (fluence rate FR, temperature T and fractional damage FD) as 3D arrays. Therefore, remember to use element-wise arithmetic operators (.* instead of *, ./ instead of / etc.)',mP_fHtrim(iM).name);
      end
      try
        if G.FRdependent(iM,2) || G.optTdependent(iM,2) || G.optFDdependent(iM,2)
          musvals = single(mP_fHtrim(iM).mus(L,FRs,Ts,FDs));
        else
          musvals = single(mP_fHtrim(iM).mus(L,0,0,0));
        end
      catch; error('Error: The mus function of medium %s is throwing an error. If you have specified the function with four input arguments, it must be able to accept the first (wavelength in nm) as a scalar and the next three (fluence rate FR, temperature T and fractional damage FD) as 3D arrays. Therefore, remember to use element-wise arithmetic operators (.* instead of *, ./ instead of / etc.)',mP_fHtrim(iM).name);
      end
      try
        if G.FRdependent(iM,3) || G.optTdependent(iM,3) || G.optFDdependent(iM,3)
          gvals = single(mP_fHtrim(iM).g(L,FRs,Ts,FDs));
        else
          gvals = single(mP_fHtrim(iM).g(L,0,0,0));
        end
      catch; error('Error: The g function of medium %s is throwing an error. If you have specified the function with four input arguments, it must be able to accept the first (wavelength in nm) as a scalar and the next three (fluence rate FR, temperature T and fractional damage FD) as 3D arrays. Therefore, remember to use element-wise arithmetic operators (.* instead of *, ./ instead of / etc.)',mP_fHtrim(iM).name);
      end
      if any(~isfinite(muavals))
        error('The mua function of medium %s returns NaN or Inf values.',mP_fHtrim(iM).name);
      elseif any(~isfinite(musvals))
        error('The mus function of medium %s returns NaN or Inf values.',mP_fHtrim(iM).name);
      elseif any(~isfinite(gvals)) && ~all(~isfinite(gvals))
        error('The g function of medium %s returns NaN or Inf values in some, but not all voxels.',mP_fHtrim(iM).name);
      end
      if all(~isfinite(gvals)) % All g values are non-finite so we use custom phase function
        useg = false;
        nThetas = 200;
        thetas = pi/nThetas*((1:nThetas) - 0.5);
        for iTheta = nThetas:-1:1
          PDF(iTheta) = sin(thetas(iTheta)).*mP_fHtrim(iM).customPhaseFunc(L,thetas(iTheta)); % The probability density function is the phase function multiplied by the differential solid angle, sin(theta)
        end
        if ~isreal(PDF) || any(~isfinite(PDF)) || any(PDF < 0)
          error('Error: The custom phase function of %s did not return real, finite, non-negative numbers.',mP_fHtrim(iM).name);
        end
        CDF = [0 cumsum(PDF)/sum(PDF)].'; % Cumulative distribution function
        matchFound = false;
        for iCDF=1:size(uniqueCDFs,2)
          if isequal(CDF,uniqueCDFs(:,iCDF))
            CDFidx = iCDF - 1; % 0-based indexing for later use in mex function
            matchFound = true;
            break;
          end
        end
        if ~matchFound
          uniqueCDFs = [uniqueCDFs CDF];
          CDFidx = size(uniqueCDFs,2) - 1; % 0-based indexing for later use in mex function
        end
      else
        useg = true;
        CDFidx = NaN;
      end

      n = mP_fHtrim(iM).n(L);
      if ~isreal(n) || any(isnan(n)) || any(n < 1)
        error('Error: The refractive index function of %s did not return real, non-NaN numbers greater than or equal to 1.',mP_fHtrim(iM).name);
      end

      nSM = mP_fHtrim(iM).nBins; % Number of sub-media to split the iM'th medium into
      if DC.minMediumIdxToConsider == uniqueMedia(iM)
        DC.minSubmediaIdx = iSM;
      end
      if DC.maxMediumIdxToConsider == uniqueMedia(iM)
        DC.maxSubmediaIdx = iSM + nSM - 1;
      end
      if isscalar(muavals) && isscalar(musvals) && isscalar(gvals)
        % No downsampling necessary
        if iL == 1
          mP.name  (iSM:iSM+nSM-1,iL) = {mP_fHtrim(iM).name};
          M(idxs3D) = iSM;
        end
        mP.mua   (iSM:iSM+nSM-1,iL) = muavals;
        mP.mus   (iSM:iSM+nSM-1,iL) = musvals;
        if useg
          mP.g   (iSM:iSM+nSM-1,iL) = gvals;
        end
        mP.n     (iSM:iSM+nSM-1,iL) = double(n);
        mP.CDFidx(iSM:iSM+nSM-1,iL) = CDFidx;
      else
        if iL == 1 % As a quick-and-dirty solution to broadband FRTFD dependent simulation binning, we just determine the binning with the first wavelength
          paramVals = NaN(3,sum(idxs3D),'single'); % Number of columns is number of voxels with this medium, and rows are (mua, mus, g)
          paramVals(1,:) = muavals;
          paramVals(2,:) = musvals;
          if useg
            paramVals(3,:) = gvals;
          else
            paramVals(3,:) = 0;
          end
          [downsampledParamVals,binidx] = getDownsampledParamVals(paramVals,nSM); % paramVals and downsampledParamVals are single precision and binidx is uint8
          [~,sortidx] = sortrows(downsampledParamVals.');
          downsampledParamVals = downsampledParamVals(:,sortidx.'); % Sort the downsampled parameter values
          I_rev = zeros(nSM,1);
          I_rev(sortidx) = 1:nSM;
          binidx = I_rev(binidx); % binidx is rearranged according to the same sort

          M(idxs3D) = iSM + binidx - 1;
          mP.name  (iSM:iSM+nSM-1  ) = {deal(mP_fHtrim(iM).name)};
          mP.mua   (iSM:iSM+nSM-1,1) = double(downsampledParamVals(1,:));
          mP.mus   (iSM:iSM+nSM-1,1) = double(downsampledParamVals(2,:));
          if useg
            mP.g   (iSM:iSM+nSM-1,1) = double(downsampledParamVals(3,:));
          end
          mP.n     (iSM:iSM+nSM-1,1) = double(n);
          mP.CDFidx(iSM:iSM+nSM-1,1) = double(CDFidx);
        else % Aside from the first wavelength, we just set the mua, mus and g for all the subsequent wavelengths to the means of the new values in the previously binned voxels
          for k=0:nSM-1
            if any(M(:) == iSM + k)
              if isscalar(muavals)
                mP.mua(iSM+k,iL) = double(muavals);
              else
                mP.mua(iSM+k,iL) = mean(double(muavals(binidx == k+1))); % Convert to double before taking the mean to reduce the accumulation of numeric errors in the addition inside the mean function
              end
              if isscalar(musvals)
                mP.mus(iSM+k,iL) = double(musvals);
              else
                mP.mus(iSM+k,iL) = mean(double(musvals(binidx == k+1))); % Convert to double before taking the mean to reduce the accumulation of numeric errors in the addition inside the mean function
              end
              if isscalar(gvals)
                mP.g(iSM+k,iL) = double(gvals);
              else
                mP.g(iSM+k,iL) = mean(double(gvals(binidx == k+1))); % Convert to double before taking the mean to reduce the accumulation of numeric errors in the addition inside the mean function
              end
            else
              mP.mua(iSM+k,iL)    = mP.mua(iSM,iL);
              mP.mus(iSM+k,iL)    = mP.mus(iSM,iL);
              mP.g  (iSM+k,iL)    = mP.g  (iSM,iL);
            end
            mP.n(iSM+k,iL)      = double(n);
            mP.CDFidx(iSM+k,iL) = double(CDFidx);
          end
        end
      end
    end
    iSM = iSM + nSM;
  end

  %% If simulating with matched interfaces, we just set all refractive indices to 1
  if matchedInterfaces
    mP.n(:) = 1;
  end

  %% Throw an error if a variable doesn't conform to its required interval
  if any(~isfinite(mP.mua(:))) || any(mP.mua(:) <= 0) || ~isreal(mP.mua)
    error('Error: A medium has invalid mua. Must be finite, real and positive.');
  elseif any(~isfinite(mP.mus(:))) || any(mP.mus(:) <= 0) || ~isreal(mP.mus)
    error('Error: A medium has invalid mus. Must be finite, real and positive.');
  elseif ~all(isnan(mP.g(:)) | (isfinite(mP.g(:)) & abs(mP.g(:)) <= 1)) || ~isreal(mP.g)
    error('Error: A medium has invalid g. Must be finite, real and between -1 and 1, or NaN.');
  elseif any(isnan(mP.n(:))) || any(mP.n(:) < 1) || ~isreal(mP.n)
    error('Error: A medium has invalid n. Must not be NaN, must be real and at least 1.');
  end

  %% Extract the refractive indices if not assuming matched interfaces, otherwise assume all 1's
  if ~matchedInterfaces
    [RI_unq , ~ , MtoRImap] = unique(mP.n,'rows');
    Gx = NaN(G.nx,G.ny,G.nz);
    Gy = NaN(G.nx,G.ny,G.nz);
    Gz = NaN(G.nx,G.ny,G.nz);
    for iRI_unq = 1:numel(RI_unq)
      n_mat = MtoRImap(M) == iRI_unq;
      [Gy_Sobel, Gx_Sobel, Gz_Sobel] = imgradientxyz(n_mat);
      if smoothingLengthScale > 0
        G_cell = smoothn({Gx_Sobel/G.dx,Gy_Sobel/G.dy,Gz_Sobel/G.dz},max(1,round(smoothingLengthScale/max([G.dx,G.dy,G.dz]))),struct('Spacing',[G.dx G.dy G.dz]));
      else
        G_cell = {Gx_Sobel/G.dx,Gy_Sobel/G.dy,Gz_Sobel/G.dz};
      end
      Gx(n_mat) = G_cell{1}(n_mat);
      Gy(n_mat) = G_cell{2}(n_mat);
      Gz(n_mat) = G_cell{3}(n_mat);
    end
    [phi,elevation,~] = cart2sph(Gx,Gy,Gz);
    theta = pi/2 - elevation;
    interfaceNormals = single([theta(:).' ; phi(:).']);
  else
    interfaceNormals = single(NaN);
  end

  if simType == 1
    model.MC.mediaProperties = mP;
    model.MC.interfaceNormals = interfaceNormals;
    model.MC.M = M;
    model.MC.CDFs = uniqueCDFs;
    model.MC.depositionCriteria = DC;
  else
    model.FMC.mediaProperties = mP;
    model.FMC.interfaceNormals = interfaceNormals;
    model.FMC.M = M;
    model.FMC.CDFs = uniqueCDFs;
    model.FMC.depositionCriteria = DC;
  end
end