function model = runMonteCarlo(model,varargin)
  %   Requires
  %       MCmatlab.mex (architecture specific)
  %
  %   See also plot, runMonteCarloFluorescence, simulateHeatDistribution

  %%%%%
  %   Copyright 2017, 2018 by Dominik Marti and Anders K. Hansen
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

  % checking on macOS whether the mex-file is quarantined and clearing it
  if ismac
    [~,cmdout] = system('xattr -l ./+MCmatlab/@model/private/MCmatlab.mexmaci64');
    if contains(cmdout,'com.apple.quarantine'); system('xattr -d com.apple.quarantine ./+MCmatlab/@model/private/MCmatlab.mexmaci64'); end
    [~,cmdout] = system('xattr -l ./+MCmatlab/@model/private/getDownsampledParamVals.mexmaci64');
    if contains(cmdout,'com.apple.quarantine'); system('xattr -d com.apple.quarantine ./+MCmatlab/@model/private/getDownsampledParamVals.mexmaci64'); end
  end

  forceSingleThreaded = isunix && ~ismac && verLessThan('matlab','9.9');
  if forceSingleThreaded % Multithreading doesn't work properly on Linux older than R2020b for some reason
    warning('Linux MCmatlab multithreading is disabled on MATLAB versions older than R2020b. Update your MATLAB to utilize multithreading.');
  end

  model.G = model.G.updateGeometry;
  G = model.G;

  if any(strcmp(varargin,'fluorescence'))
    simType = 2;
    useGPU = model.FMC.useGPU;
  else
    simType = 1;
    useGPU = model.MC.useGPU;
    % Calculate spectrum from the spectrum function handle:
    for iWavelength = numel(model.MC.wavelength):-1:1
      model.MC.spectrum(iWavelength) = model.MC.spectrumFunc(model.MC.wavelength(iWavelength));
    end
    model.MC.spectrum = model.MC.spectrum/sum(model.MC.spectrum); % Normalize
    model.MC.sourceDistribution = model.MC.sourceDistribution/(model.G.dx*model.G.dy*model.G.dz*sum(model.MC.sourceDistribution(:))); % Normalize
  end

  checkMCinputFields(model,simType);

  %% Get initial temperature, fractional damage and fluence rate
  T = NaN([G.nx G.ny G.nz],'single');
  T(:) = model.HS.T; % model.HS.T may be scalar or 3D

  FD = NaN([G.nx G.ny G.nz],'single');
  FD(:) = 1 - exp(-model.HS.Omega); % Fractional damage of molecules/cells. model.HS.Omega may be scalar 0 or 3D

  FR = NaN([G.nx G.ny G.nz],'single');
  FR(:) = model.MC.P*sum(model.MC.NFR,4); % All excitation wavelengths are summed to get the total fluence rate. model.MC.NFR may be scalar 0 or 3D or 4D

  %% Call Monte Carlo C code (MEX file) to get fluence rate distribution, possibly inside an iterative loop in the case of fluence rate dependent excitation light
  iterate = simType == 1 && model.MC.FRdepIterations > 0;

  if iterate
    nIterations = model.MC.FRdepIterations;
    if isnan(model.MC.nPhotonsRequested)
      t = model.MC.simulationTimeRequested;
      t_array = t*(1/2).^(nIterations-1:-1:0);
    else
      nPhotons = model.MC.nPhotonsRequested;
      nPhotons_array = nPhotons*(1/2).^(nIterations-1:-1:0);
    end
    for i=1:nIterations
      if isnan(model.MC.nPhotonsRequested)
        model.MC.simulationTimeRequested = t_array(i);
      else
        model.MC.nPhotonsRequested = nPhotons_array(i);
      end
      model = getOpticalMediaProperties(model,simType); % Also performs splitting of mediaProperties and M_raw, if necessary
      if useGPU
        model = MCmatlab_CUDA(model,simType);
      else
        if forceSingleThreaded % Multithreading doesn't work properly on Linux older than R2020b for some reason
          model = MCmatlab_singlethreaded(model,simType);
        else
          model = MCmatlab(model,simType);
        end
      end
      if i<nIterations; model.MC.FR = model.MC.P*model.MC.NFR; end
    end
    if isnan(model.MC.nPhotonsRequested)
      model.MC.simulationTimeRequested = t;
    else
      model.MC.nPhotonsRequested = nPhotons;
    end
  else
    model = getOpticalMediaProperties(model,simType); % Also performs splitting of mediaProperties and M_raw, if necessary
    if simType == 2 % Since we're not iterating, we might be simulating fluorescence
      %% Calculate 3D fluorescence source distribution
      Ls_E = model.MC.wavelength; % Excitation wavelengths
      nL_E = numel(Ls_E); % Number of excitation wavelengths
      Ls_F = model.FMC.wavelength; % Fluorescence wavelengths
      nL_F = numel(Ls_F); % Number of fluorescence wavelengths
      mP_fH = G.mediaPropertiesFunc(G.mediaPropParams);
      nM = numel(mP_fH);
      Fphotons = zeros(size(G.M_raw)); % Proportional to how many fluorescence photons are emitted, as function of xyz, summed over all wavelengths
      emSp = NaN([size(G.M_raw),nL_F]); % The emission spectra as function of wavelength. These start out non-normalized but will be normalized later. Because each emission spectrum may depend on FR, T and FD, each voxel may have its own emission spectrum
      NA = model.MC.NA;
      for iM = 1:nM
        idxs = find(G.M_raw == iM);
        FRs = FR(idxs);
        Ts = T(idxs);
        FDs = FD(idxs);
        Fphotons_thisMedium = zeros(size(idxs));
        for iL_E = 1:nL_E % Loop over excitation wavelengths
          try
            if G.FRdependent(iM,4) || G.optTdependent(iM,4) || G.optFDdependent(iM,4)
              QYs = mP_fH(iM).QY(Ls_E(iL_E),FRs,Ts,FDs);
            else
              QYs = mP_fH(iM).QY(Ls_E(iL_E),0,0,0);
            end
          catch; error('Error: The QY function of medium %s is throwing an error. If you have specified the function with four input arguments, it must be able to accept the first (wavelength in nm) as a scalar and the next three (fluence rate FR, temperature T and fractional damage FD) as 3D arrays. Therefore, remember to use element-wise arithmetic operators (.* instead of *, ./ instead of / etc.)',mP_fH(iM).name);
          end
          Fphotons_thisMedium = Fphotons_thisMedium + NA(idxs + (iL_E-1)*G.nx*G.ny*G.nz).*QYs.*Ls_E(iL_E);
        end
        Fphotons(idxs) = Fphotons_thisMedium;
        for iL_F = 1:nL_F % Loop over fluorescence wavelengths
          try 
            if G.FRdependent(iM,5) || G.optTdependent(iM,5) || G.optFDdependent(iM,5)
              emSp(idxs + (iL_F-1)*G.nx*G.ny*G.nz) = mP_fH(iM).ES(Ls_F(iL_F),FRs,Ts,FDs);
            else
              emSp(idxs + (iL_F-1)*G.nx*G.ny*G.nz) = mP_fH(iM).ES(Ls_F(iL_F),0,0,0);
            end
          catch; error('Error: The ES function of medium %s is throwing an error. If you have specified the function with four input arguments, it must be able to accept the first (wavelength in nm) as a scalar and the next three (fluence rate FR, temperature T and fractional damage FD) as 3D arrays. Therefore, remember to use element-wise arithmetic operators (.* instead of *, ./ instead of / etc.)',mP_fH(iM).name);
          end
        end
      end
      if any(~isfinite(Fphotons(:))) || ~isreal(Fphotons) || any(Fphotons(:) < 0)
        error('Error: A quantum yield function is returning values that are not finite, real and non-negative.');
      end
      if any(~isfinite(emSp(:))) || ~isreal(emSp) || any(emSp(:) < 0)
        error('Error: An emission spectrum function is returning values that are not finite, real and non-negative.');
      end
      
      emSp = emSp./repmat(sum(emSp,4),1,1,1,nL_F); % Normalization
      model.FMC.sourceDistribution = repmat(Fphotons,1,1,1,nL_F).*emSp./repmat(permute(Ls_F(:),[2 3 4 1]),G.nx,G.ny,G.nz,1); % The fluorescence source distribution is each voxel's emission spectrum normalized to a sum equal to the corresponding value in Fsources3D

      if max(model.FMC.sourceDistribution(:)) == 0; error('Error: No fluorescence emitters'); end
    end
    if useGPU
      model = MCmatlab_CUDA(model,simType);
    else
      if forceSingleThreaded % Multithreading doesn't work properly on Linux older than R2020b for some reason
        model = MCmatlab_singlethreaded(model,simType);
      else
        model = MCmatlab(model,simType);
      end
    end
  end

  if simType == 2
    % Add positions of the centers of the pixels in the light collector
    % image and the time array
    if model.FMC.useLightCollector && model.FMC.LC.res > 1
      model.FMC.LC.X = linspace(model.FMC.LC.fieldSize*(1/model.FMC.LC.res-1),model.FMC.LC.fieldSize*(1-1/model.FMC.LC.res),model.FMC.LC.res)/2;
      model.FMC.LC.Y = model.FMC.LC.X;
    end
    if model.FMC.useLightCollector && model.FMC.LC.nTimeBins > 0
      if model.FMC.matchedInterfaces
        warning('Time tagging is on, but matchedInterfaces is true, which sets all refractive indices to 1.');
      end
      model.FMC.LC.t = (-1/2:(model.FMC.LC.nTimeBins+1/2))*(model.FMC.LC.tEnd-model.FMC.LC.tStart)/model.FMC.LC.nTimeBins + model.FMC.LC.tStart;
    end

    % Add angles of the centers of the far field pixels
    if model.FMC.farFieldRes
      model.FMC.farFieldTheta = linspace(pi/model.FMC.farFieldRes/2,pi-pi/model.FMC.farFieldRes/2,model.FMC.farFieldRes);
      model.FMC.farFieldPhi   = linspace(-pi+(2*pi)/model.FMC.farFieldRes/2,pi-(2*pi)/model.FMC.farFieldRes/2,model.FMC.farFieldRes);
    end
  else
    % Add positions of the centers of the pixels in the light collector
    % image and the time array
    if model.MC.useLightCollector && model.MC.LC.res > 1
      model.MC.LC.X = linspace(model.MC.LC.fieldSize*(1/model.MC.LC.res-1),model.MC.LC.fieldSize*(1-1/model.MC.LC.res),model.MC.LC.res)/2;
      model.MC.LC.Y = model.MC.LC.X;
    end
    if model.MC.useLightCollector && model.MC.LC.nTimeBins > 0
      if model.MC.matchedInterfaces
        warning('Time tagging is on, but matchedInterfaces is true, which sets all refractive indices to 1.');
      end
      model.MC.LC.t = (-1/2:(model.MC.LC.nTimeBins+1/2))*(model.MC.LC.tEnd-model.MC.LC.tStart)/model.MC.LC.nTimeBins + model.MC.LC.tStart;
    end

    % Add angles of the centers of the far field pixels
    if model.MC.farFieldRes
      model.MC.farFieldTheta = linspace(pi/model.MC.farFieldRes/2,pi-pi/model.MC.farFieldRes/2,model.MC.farFieldRes);
      model.MC.farFieldPhi   = linspace(-pi+(2*pi)/model.MC.farFieldRes/2,pi-(2*pi)/model.MC.farFieldRes/2,model.MC.farFieldRes);
    end
  end
end

function checkMCinputFields(model,simType)
  G = model.G;

  if simType == 2
    MCorFMC = model.FMC;
  else
    MCorFMC = model.MC;
  end

  %% Check for GPU compatibility if needed
  if (simType == 1 && model.MC.useGPU) || (simType == 2 && model.FMC.useGPU)
    if ~ispc
      error('GPU computing is only supported on Windows');
    end
    v = ver;
    if ~any(strcmp({v.Name},'Parallel Computing Toolbox'))
      error('You must have the Parallel Computing Toolbox installed to use GPU acceleration');
    end
    try
      GPUDev = gpuDevice;
      if(str2double(GPUDev.ComputeCapability) < 3)
        error('Your GPU is too old (CUDA compute capability < 3.0)')
      end
    catch
      error('No supported NVIDIA GPU found, or its driver is too old');
    end
  end

  %% Check other fields
  if MCorFMC.depositionCriteria.onlyCollected && ~MCorFMC.depositionCriteria.evaluateOnlyAtEndOfLife
    error('Error: depositionCriteria.onlyCollected = true requires depositionCriteria.evaluateOnlyAtEndOfLife = true.');
  end
  if MCorFMC.nExamplePaths ~= 0 && ~isscalar(MCorFMC.wavelength)
    error('Error: Example paths are not supported when simulating broadband illumination.');
  end
  if isnan(MCorFMC.wavelength)
    error('Error: No wavelength defined');
  end
  if MCorFMC.depositionCriteria.onlyCollected && ~MCorFMC.useLightCollector
    error('Error: depositionCriteria.onlyCollected is true, but no light collector is defined');
  end
  if MCorFMC.farFieldRes && MCorFMC.boundaryType == 0
    error('Error: If boundaryType == 0, no photons can escape to be registered in the far field. Set farFieldRes to zero or change boundaryType.');
  end

  if simType == 2
    if isscalar(model.MC.NFR)
      error('Error: normalizedFluenceRate matrix not calculated for excitation light');
    end
    DC = model.MC.depositionCriteria;
    if DC.minScatterings ~= 0 || ~isinf(DC.maxScatterings) || ...
       DC.minRefractions ~= 0 || ~isinf(DC.maxRefractions) || ...
       DC.minReflections ~= 0 || ~isinf(DC.maxReflections) || ...
       DC.minInterfaceTransitions ~= 0 || ~isinf(DC.maxInterfaceTransitions) || ...
       DC.onlyCollected
      error('The excitation step must be calculated without restrictive deposition criteria in order to calculate fluorescence.');
    end
  else
    if ~isscalar(model.MC.sourceDistribution) % Distributed source
      size4D = size(model.MC.sourceDistribution);
      size4D(4) = size(model.MC.sourceDistribution,4);
      if ~isequal(size4D, [model.G.nx, model.G.ny, model.G.nz, numel(model.MC.wavelength)])
        error('Error: model.MC.sourceDistribution must be of size [nx, ny, nz, nLambda]');
      end
      if ~isreal(model.MC.sourceDistribution) || ~all(isfinite(model.MC.sourceDistribution(:))) || ~all(model.MC.sourceDistribution(:) >= 0)
        error('Error: model.MC.sourceDistribution must be real, finite and non-negative');
      end
    else % Not distributed source
      if isnan(MCorFMC.LS.sourceType)
        error('Error: No sourceType defined');
      end
      if MCorFMC.LS.sourceType == 4 && (isnan(MCorFMC.LS.FPID.radialDistr(1)) || isnan(MCorFMC.LS.AID.radialDistr(1)))
        error('Error: lightSource.focalPlaneIntensityDistribution.radialDistr and lightSource.angularIntensityDistribution.radialDistr must both be specified when sourceType is 4');
      end
      if MCorFMC.LS.sourceType == 3
        if isnan(MCorFMC.LS.FPID.radialWidth) || isnan(MCorFMC.LS.AID.radialWidth)
          error('Error: lightSource.focalPlaneIntensityDistribution.radialWidth and lightSource.angularIntensityDistribution.radialWidth must both be specified when sourceType is 3');
        end
      end
      if MCorFMC.LS.sourceType == 4
        if isnan(MCorFMC.LS.FPID.radialWidth)
          error('Error: lightSource.focalPlaneIntensityDistribution.radialWidth must be specified when sourceType is 4');
        end      
        if ~isequal(MCorFMC.LS.AID.radialDistr, 2) && isnan(MCorFMC.LS.FPID.radialWidth) % If non-Lambertian source
          error('Error: lightSource.angularIntensityDistribution.radialWidth must be specified when sourceType is 4 and the source is non-Lambertian');
        end
      end
      if MCorFMC.LS.sourceType == 5
        if isnan(MCorFMC.LS.FPID.XDistr(1))
          error('Error: lightSource.focalPlaneIntensityDistribution.XDistr must be specified when sourceType is 5');
        end
        if isnan(MCorFMC.LS.FPID.XWidth)
          error('Error: lightSource.focalPlaneIntensityDistribution.Xwidth must be specified when sourceType is 5');
        end
        if isnan(MCorFMC.LS.FPID.YDistr(1))
          error('Error: lightSource.focalPlaneIntensityDistribution.YDistr must be specified when sourceType is 5');
        end
        if isnan(MCorFMC.LS.FPID.YWidth)
          error('Error: lightSource.focalPlaneIntensityDistribution.YWidth must be specified when sourceType is 5');
        end
        if isnan(MCorFMC.LS.AID.XDistr(1))
          error('Error: lightSource.angularIntensityDistribution.XDistr must be specified when sourceType is 5');
        end
        if MCorFMC.LS.AID.XDistr(1) ~= 2 && isnan(MCorFMC.LS.AID.XWidth)
          error('Error: lightSource.angularIntensityDistribution.Xwidth must be specified when sourceType is 5 and the distribution is non-Lambertian');
        end
        if isnan(MCorFMC.LS.AID.YDistr(1))
          error('Error: lightSource.angularIntensityDistribution.YDistr must be specified when sourceType is 5');
        end
        if MCorFMC.LS.AID.YDistr(1) ~= 2 && isnan(MCorFMC.LS.AID.YWidth)
          error('Error: lightSource.angularIntensityDistribution.YWidth must be specified when sourceType is 5 and the distribution is non-Lambertian');
        end
        if xor(MCorFMC.LS.AID.XDistr == 2,MCorFMC.LS.AID.YDistr == 2)
          error('Error: lightSource.angularIntensityDistribution.XDistr and lightSource.angularIntensityDistribution.YDistr must either both be set to cosine (Lambertian), or neither');
        end
      end
      if size(MCorFMC.LS.FPID.radialDistr,1) > 1 || size(MCorFMC.LS.FPID.XDistr,1) > 1 || size(MCorFMC.LS.FPID.YDistr,1) > 1 ||...
          size(MCorFMC.LS.AID.radialDistr,1) > 1 || size(MCorFMC.LS.AID.XDistr,1) > 1 || size(MCorFMC.LS.AID.YDistr,1) > 1
        error('Error: Beam distribution functions must be scalars or row-vectors, not column-vectors');
      end
      if (numel(MCorFMC.LS.FPID.radialDistr) > 1 && numel(MCorFMC.LS.FPID.radialDistr) < 1000) ||...
          (numel(MCorFMC.LS.FPID.XDistr)      > 1 && numel(MCorFMC.LS.FPID.XDistr)      < 1000) ||...
          (numel(MCorFMC.LS.FPID.YDistr)      > 1 && numel(MCorFMC.LS.FPID.YDistr)      < 1000) ||...
          (numel(MCorFMC.LS.AID.radialDistr) > 1 && numel(MCorFMC.LS.AID.radialDistr) < 1000) ||...
          (numel(MCorFMC.LS.AID.XDistr)      > 1 && numel(MCorFMC.LS.AID.XDistr)      < 1000) ||...
          (numel(MCorFMC.LS.AID.YDistr)      > 1 && numel(MCorFMC.LS.AID.YDistr)      < 1000)
        error('Error: Beam definition distributions must have 1000 elements (or more, if necessary)');
      end
    end
  end

  if MCorFMC.useLightCollector
    %% Check to ensure that the light collector is not inside the cuboid and set res to 1 if using fiber
    if MCorFMC.boundaryType == 0
      error('Error: If boundaryType == 0, no photons can escape to be registered on the light collector. Disable light collector or change boundaryType.');
    end
    if isfinite(MCorFMC.LC.f)
      xLCC = MCorFMC.LC.x - MCorFMC.LC.f*sin(MCorFMC.LC.theta)*cos(MCorFMC.LC.phi); % x position of Light Collector Center
      yLCC = MCorFMC.LC.y - MCorFMC.LC.f*sin(MCorFMC.LC.theta)*sin(MCorFMC.LC.phi); % y position
      zLCC = MCorFMC.LC.z - MCorFMC.LC.f*cos(MCorFMC.LC.theta);                     % z position
    else
      xLCC = MCorFMC.LC.x;
      yLCC = MCorFMC.LC.y;
      zLCC = MCorFMC.LC.z;
      if MCorFMC.LC.res ~= 1
        error('Error: lightCollector.res must be 1 when lightCollector.f is Inf');
      end
    end

    if (abs(xLCC)               < G.nx*G.dx/2 && ...
        abs(yLCC)               < G.ny*G.dy/2 && ...
        abs(zLCC - G.nz*G.dz/2) < G.nz*G.dz/2)
      error('Error: Light collector center (%.4f,%.4f,%.4f) is inside cuboid',xLCC,yLCC,zLCC);
    end
  end
end