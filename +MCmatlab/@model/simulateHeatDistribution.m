function model = simulateHeatDistribution(model)
  %   Simulates temperature evolution due to absorption of light and
  %   diffusion of heat through the cuboid based on output of runMonteCarlo.m.
  %   Also calculates Arrhenius-based thermal damage.
  %

  % checking on macOS whether the mex-file is quarantined and clearing it
  if ismac
    [~,cmdout] = system('xattr -l ./+MCmatlab/@model/private/finiteElementHeatPropagator.mexmaci64');
    if contains(cmdout,'com.apple.quarantine'); system('xattr -d com.apple.quarantine ./+MCmatlab/@model/private/finiteElementHeatPropagator.mexmaci64'); end
  end

  forceSingleThreaded = isunix && ~ismac && verLessThan('matlab','9.9');
  if forceSingleThreaded % Multithreading doesn't work properly on Linux older than R2020b for some reason
    warning('Linux MCmatlab multithreading is disabled on MATLAB versions older than R2020b. Update your MATLAB to utilize multithreading.');
  end

  if isscalar(model.MC.NFR)
    error('Error: MC.NFR was not calculated');
  end

  model.G = model.G.updateGeometry;
  G = model.G;

  T = NaN(size(G.M_raw),'single');
  T(:) = model.HS.T;
  model.HS.T = T;

  firstHSrun = isnan(model.HS.maxMediaTemps(1));

  model = getThermalMediaProperties(model);
  mP = model.HS.mediaProperties; % Split media

  mP_fH = model.G.mediaPropertiesFunc(G.mediaPropParams); % Unsplit media
  mP_fHtrim = mP_fH(unique(G.M_raw));
  nM = numel(mP_fHtrim);
  names = {mP_fHtrim.name};

  if firstHSrun
    model.HS.maxMediaTemps = ones(nM,1);
    for i=1:nM
      if any(G.M_raw == i)
        model.HS.maxMediaTemps(i) = max(model.HS.T(G.M_raw == i));
      else
        model.HS.maxMediaTemps(i) = -Inf;
      end
    end
  end

  %% Prepare Arrhenius thermal damage parameters
  A = mP.A;
  E = mP.E;
  if any(A) % If non-zero Arrhenius data exists, prepare to calculate thermal damage.
    if firstHSrun
      model.HS.Omega = zeros(size(G.M_raw),'single');
    end
  else
    model.HS.Omega = single(NaN);
  end

  %% Calculate amount and duration of time steps in each phase
  % If either pulse on-duration or off-duration is 0, it doesn't make sense to talk of multiple pulses
  if model.HS.nPulses ~= 1 && (model.HS.durationOn == 0 || model.HS.durationOff == 0)
    error('Error: If either on-duration or off-duration is 0, nPulses must be 1');
  end

  if model.HS.silentMode && ~any(G.optTdependent(:)) && ~any(G.thmTdependent(:)) && ~any(G.optFDdependent(:)) && ~any(G.thmFDdependent(:))
    model.HS.nUpdates = 1;
  end

  VHC = single(mP.VHC); % Volumetric heat capacity array. Element i is the VHC of medium i in the split mediaProperties list.
  TC = single(mP.TC); % Thermal conductivity array. Element i is the TC of medium i in the split mediaProperties list.
  if model.HS.largeTimeSteps
    dtmax = calcdtmax(model.HS.M,TC,VHC,G.dx,G.dy,G.dz)*15; % Highest allowable time step duration for low precision (factor 15 is found to be the highest value for which there are no artifacts visible in the example 3 test case)
  else
    dtmax = calcdtmax(model.HS.M,TC,VHC,G.dx,G.dy,G.dz)*1.5; % Highest allowable time step duration for normal precision
  end

  if model.HS.durationOn ~= 0
    nUpdatesOn = max(1,round(model.HS.nUpdates*model.HS.durationOn/(model.HS.durationOn+model.HS.durationOff)));
    nTsPerUpdateOn = ceil(model.HS.durationOn/nUpdatesOn/dtmax); % Number of time steps per update in illumination phase
    dtOn = model.HS.durationOn/nUpdatesOn/nTsPerUpdateOn; % Time step size during illumination
  else
    nUpdatesOn = 0;
    nTsPerUpdateOn = 1;
    dtOn = 1;
  end
  nDigitsOn = 1+floor(log10(nUpdatesOn));

  if model.HS.durationOff ~= 0
    nUpdatesOff = max(1,model.HS.nUpdates - nUpdatesOn);
    nTsPerUpdateOff = ceil(model.HS.durationOff/nUpdatesOff/dtmax); % Number of time steps per update in diffusion phase
    dtOff = model.HS.durationOff/nUpdatesOff/nTsPerUpdateOff; % Time step size during diffusion
  else
    nUpdatesOff = 0;
    nTsPerUpdateOff = 1;
    dtOff = 1;
  end
  nDigitsOff = 1+floor(log10(nUpdatesOff));

  if model.HS.durationEnd ~= 0
    nUpdatesEnd = max(1,round(model.HS.nUpdates*model.HS.durationEnd/(model.HS.durationOn+model.HS.durationOff)));
    nTsPerUpdateEnd = ceil(model.HS.durationEnd/nUpdatesEnd/dtmax); % Number of time steps per update in end relaxation phase
    dtEnd = model.HS.durationEnd/nUpdatesEnd/nTsPerUpdateEnd; % Time step size during end relaxation phase
  else
    nUpdatesEnd = 0;
    nTsPerUpdateEnd = 1;
    dtEnd = 1;
  end
  nDigitsEnd = 1+floor(log10(nUpdatesEnd));

  % If a previous heat simulation was performed, we need to offset the time
  if ~firstHSrun
    t_offset = model.HS.sensorsTimeVector(end);
  else
    t_offset = 0;
  end

  % Array of times associated with the updates
  updatesTimeVector = t_offset + [0 , ((model.HS.durationOn+model.HS.durationOff)*repelem(0:(model.HS.nPulses-1),nUpdatesOn                + nUpdatesOff                ) + repmat([(1:nUpdatesOn)*nTsPerUpdateOn*dtOn , (model.HS.durationOn + ((1:nUpdatesOff)*nTsPerUpdateOff*dtOff))],1,model.HS.nPulses)) , (model.HS.durationOn + model.HS.durationOff)*model.HS.nPulses + (1:nUpdatesEnd)*nTsPerUpdateEnd*dtEnd];

  % Array of times associated with the steps
  sensorsTimeVector = t_offset + [0 , ((model.HS.durationOn+model.HS.durationOff)*repelem(0:(model.HS.nPulses-1),nUpdatesOn*nTsPerUpdateOn + nUpdatesOff*nTsPerUpdateOff) + repmat([(1:nUpdatesOn*nTsPerUpdateOn)*dtOn , (model.HS.durationOn + ((1:nUpdatesOff*nTsPerUpdateOff)*dtOff))],1,model.HS.nPulses)) , (model.HS.durationOn + model.HS.durationOff)*model.HS.nPulses + (1:nUpdatesEnd*nTsPerUpdateEnd)*dtEnd];
  if isnan(model.HS.sensorsTimeVector(1))
    model.HS.sensorsTimeVector = sensorsTimeVector;
  else
    model.HS.sensorsTimeVector = [model.HS.sensorsTimeVector(1:end-1) , sensorsTimeVector];
  end

  %% Calculate proportionality between voxel-to-voxel temperature difference DeltaT and time step temperature change dT, ...
  %  as well as temperature change due to absorbed heat per time step
  [dTdtperdeltaT , dTdt_abs] = calcdTdtArrays(model);

  %% Convert the temperature sensor positions to interpolation corner indices and weights used in the MEX function
  numTemperatureSensors = size(model.HS.tempSensorPositions,1);
  if numTemperatureSensors
    if any(abs(model.HS.tempSensorPositions(:,1))               > G.dx*G.nx/2) || ...
        any(abs(model.HS.tempSensorPositions(:,2))               > G.dy*G.ny/2) || ...
        any(abs(model.HS.tempSensorPositions(:,3)) - G.dz*G.nz/2 > G.dz*G.nz/2)
      error('Error: A temperature sensor is outside the cuboid.');
    end

    [tempSensorCornerIdxs , tempSensorInterpWeights] = convertTempSensorPos(model.HS.tempSensorPositions,G.x,G.y,G.z);
  else
    tempSensorCornerIdxs = [];
    tempSensorInterpWeights = [];
  end
  if isnan(model.HS.sensorTemps(1))
    model.HS.sensorTemps = [];
  end

  %% Output some diagnostics info and prepare the temperature evolution plot. If making a movie, put a geometry illustration into the beginning of the movie.
  if ~model.HS.silentMode
    fprintf('---------------------Heat Simulation---------------------\n');
    if model.HS.durationOn ~= 0
      fprintf('Illumination phase consists of %d steps of %0.2e s.\n',nUpdatesOn*nTsPerUpdateOn,dtOn);
    end
    if model.HS.durationOff ~= 0
      fprintf('Diffusion phase consists of %d steps of %0.2e s.\n',nUpdatesOff*nTsPerUpdateOff,dtOff);
    end
    if model.HS.durationEnd ~= 0
      fprintf('End phase consists of %d steps of %0.2e s.\n',nUpdatesEnd*nTsPerUpdateEnd,dtEnd);
    end

    movieFrameidx = 1;
    heatsimFigure = figure(21);
    if model.HS.makeMovie
      set(heatsimFigure,'WindowStyle','Normal');
      drawnow;
      heatsimFigure.Position = [40 160 1100 650]; % The size defined here will determine the movie resolution
      if firstHSrun % Make a temporary figure showing the geometry illustration to put into the beginning of the movie
        heatsimFigure = MCmatlab.NdimSliderPlot(G.M_raw,...
          'nFig',21,...
          'axisValues',{G.x,G.y,G.z},...
          'axisLabels',{'x [cm]','y [cm]','z [cm]','Geometry illustration'},...
          'indexLabels',names,...
          'linColormap',lines,...
          'axisEqual',true,...
          'reversedAxes',3,...
          'slicePositions',model.HS.slicePositions,...
          'docked',~model.HS.makeMovie);
        drawnow;
        movieFrames(movieFrameidx) = getframe(heatsimFigure);
        movieFrameidx = movieFrameidx + 1;
      end
    else
      set(heatsimFigure,'WindowStyle','Docked');
    end

    heatsimFigure = MCmatlab.NdimSliderPlot(model.HS.T,...
      'nFig',21,...
      'axisValues',{G.x,G.y,G.z},...
      'axisLabels',{'x [cm]','y [cm]','z [cm]','Temperature evolution'},...
      'linColormap',hot(256),...
      'axisEqual',true,...
      'reversedAxes',3,...
      'slicePositions',model.HS.slicePositions,...
      'docked',~model.HS.makeMovie);
    heatsimFigure.Name = 'Temperature evolution';
    h_title = title(['Temperature evolution, t = ' num2str(updatesTimeVector(1),'%#.2g') ' s']);
    caxis(model.HS.plotTempLimits); % User-defined color scale limits
    if model.HS.makeMovie && firstHSrun
      movieFrames(movieFrameidx) = getframe(heatsimFigure);
      movieFrameidx = movieFrameidx + 1;
    end
  end

  %% Put heatSim parameters into a struct
  heatSimParameters = struct('M',model.HS.M-1,'A',A,'E',E,'dTdtperdeltaT',dTdtperdeltaT,'dTdt_abs',dTdt_abs,...
    'useAllCPUs',model.HS.useAllCPUs,'heatBoundaryType',model.HS.heatBoundaryType,...
    'tempSensorCornerIdxs',tempSensorCornerIdxs,'tempSensorInterpWeights',tempSensorInterpWeights); % Contents of HS.M have to be converted from Matlab's 1-based indexing to C's 0-based indexing.

  %% Check for GPU compatibility if needed, and start timer
  if model.HS.useGPU
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
  if ~model.HS.silentMode
    if model.HS.useGPU
      GPUDevice = gpuDevice;
      fprintf("Using %s with CUDA compute capability %s\n",GPUDevice.Name,GPUDevice.ComputeCapability);
    end
    tic;
  end

  %% Simulate heat transfer
  for j=1:model.HS.nPulses
    %% Illumination phase
    if model.HS.durationOn ~= 0
      if ~model.HS.silentMode
        updateIdx = (j-1)*(nUpdatesOn+nUpdatesOff); % Index of the update
        if updateIdx > 0
          timeRemainingEstimate = toc/updateIdx*((nUpdatesOn+nUpdatesOff)*model.HS.nPulses+nUpdatesEnd - updateIdx);
          nCharsOnLine = fprintf(1,'[%02.0f:%02.0f:%02.0f] Illuminating pulse #%d... 0/%d',floor(timeRemainingEstimate/3600),floor(rem(timeRemainingEstimate/60,60)),floor(rem(timeRemainingEstimate,60)),j,nUpdatesOn);
        else
          nCharsOnLine = fprintf(1,'[00:00:00] Illuminating pulse #%d... 0/%d',j,nUpdatesOn);
        end
        drawnow;
      end

      heatSimParameters.lightsOn = true;
      heatSimParameters.steps = nTsPerUpdateOn;
      heatSimParameters.dt = dtOn;
      for i = 1:nUpdatesOn
        if model.HS.useGPU
          [model.HS.T,model.HS.Omega,newSensorTemps,newMaxMediaTemps] = finiteElementHeatPropagator_CUDA(model.HS.T,model.HS.Omega,heatSimParameters);
        else
          if forceSingleThreaded % Multithreading doesn't work properly on Linux older than R2020b for some reason
            [model.HS.T,model.HS.Omega,newSensorTemps,newMaxMediaTemps] = finiteElementHeatPropagator_singlethreaded(model.HS.T,model.HS.Omega,heatSimParameters);
          else
            [model.HS.T,model.HS.Omega,newSensorTemps,newMaxMediaTemps] = finiteElementHeatPropagator(model.HS.T,model.HS.Omega,heatSimParameters);
          end
        end

        SMidx = 1; % sub-medium index
        for iM=1:nM
          model.HS.maxMediaTemps(iM) = max([model.HS.maxMediaTemps(iM), double(newMaxMediaTemps(SMidx:(SMidx + mP_fHtrim(iM).nBins - 1)))]);
          SMidx = SMidx + mP_fHtrim(iM).nBins;
        end

        if isempty(model.HS.sensorTemps)
          model.HS.sensorTemps = newSensorTemps;
        else
          model.HS.sensorTemps = [model.HS.sensorTemps(:,1:end-1) newSensorTemps];
        end

        updateIdx = i+(j-1)*(nUpdatesOn+nUpdatesOff); % Index of the update
        if ~model.HS.silentMode
          fprintf(repmat('\b',1,nCharsOnLine));
          timeRemainingEstimate = toc/updateIdx*((nUpdatesOn+nUpdatesOff)*model.HS.nPulses+nUpdatesEnd - updateIdx);
          nCharsOnLine = fprintf(1,'[%02.0f:%02.0f:%02.0f] Illuminating pulse #%d... %d/%d',floor(timeRemainingEstimate/3600),floor(rem(timeRemainingEstimate/60,60)),floor(rem(timeRemainingEstimate,60)),j,i,nUpdatesOn);
          updateTempPlot(heatsimFigure,model.HS.T);
          h_title.String = ['Temperature evolution, t = ' num2str(updatesTimeVector(updateIdx+1),'%#.2g') ' s'];
          drawnow;

          if model.HS.makeMovie
            movieFrames(movieFrameidx) = getframe(heatsimFigure);
            movieFrameidx = movieFrameidx + 1;
          end
        end

        if rem(updateIdx,model.HS.mediaPropRecalcPeriod) == 0 % If this update is one where we have to check for recalculating the optical or thermal properties
          [model,heatSimParameters] = heatSimParamsRecalc(model,heatSimParameters); % Recalculate MC and FMC, if necessary, and get new thermal properties and a new media matrix
        end
      end
      if ~model.HS.silentMode; fprintf(1,repmat('\b',1,nCharsOnLine)); end
    end

    %% Diffusion phase
    if model.HS.durationOff ~= 0
      if ~model.HS.silentMode
        updateIdx = nUpdatesOn+(j-1)*(nUpdatesOn+nUpdatesOff); % Index of the update
        if updateIdx > 0
          timeRemainingEstimate = toc/updateIdx*((nUpdatesOn+nUpdatesOff)*model.HS.nPulses+nUpdatesEnd - updateIdx);
          nCharsOnLine = fprintf(1,'[%02.0f:%02.0f:%02.0f] Diffusing heat after pulse #%d... 0/%d',floor(timeRemainingEstimate/3600),floor(rem(timeRemainingEstimate/60,60)),floor(rem(timeRemainingEstimate,60)),j,nUpdatesOff);
        else
          nCharsOnLine = fprintf(1,'[00:00:00] Diffusing heat after pulse #%d... 0/%d',j,nUpdatesOff);
        end
        drawnow;
      end

      heatSimParameters.lightsOn = false;
      heatSimParameters.steps = nTsPerUpdateOff;
      heatSimParameters.dt = dtOff;
      for i = 1:nUpdatesOff
        if model.HS.useGPU
          [model.HS.T,model.HS.Omega,newSensorTemps,newMaxMediaTemps] = finiteElementHeatPropagator_CUDA(model.HS.T,model.HS.Omega,heatSimParameters);
        else
          if forceSingleThreaded % Multithreading doesn't work properly on Linux older than R2020b for some reason
            [model.HS.T,model.HS.Omega,newSensorTemps,newMaxMediaTemps] = finiteElementHeatPropagator_singlethreaded(model.HS.T,model.HS.Omega,heatSimParameters);
          else
            [model.HS.T,model.HS.Omega,newSensorTemps,newMaxMediaTemps] = finiteElementHeatPropagator(model.HS.T,model.HS.Omega,heatSimParameters);
          end
        end

        SMidx = 1; % sub-medium index
        for iM=1:nM
          model.HS.maxMediaTemps(iM) = max([model.HS.maxMediaTemps(iM), double(newMaxMediaTemps(SMidx:(SMidx + mP_fHtrim(iM).nBins - 1)))]);
          SMidx = SMidx + mP_fHtrim(iM).nBins;
        end

        if isempty(model.HS.sensorTemps)
          model.HS.sensorTemps = newSensorTemps;
        else
          model.HS.sensorTemps = [model.HS.sensorTemps(:,1:end-1) newSensorTemps];
        end

        updateIdx = nUpdatesOn+i+(j-1)*(nUpdatesOn+nUpdatesOff); % Index of the update
        if ~model.HS.silentMode
          fprintf(repmat('\b',1,nCharsOnLine));
          timeRemainingEstimate = toc/updateIdx*((nUpdatesOn+nUpdatesOff)*model.HS.nPulses+nUpdatesEnd - updateIdx);
          nCharsOnLine = fprintf(1,'[%02.0f:%02.0f:%02.0f] Diffusing heat after pulse #%d... %d/%d',floor(timeRemainingEstimate/3600),floor(rem(timeRemainingEstimate/60,60)),floor(rem(timeRemainingEstimate,60)),j,i,nUpdatesOff);
          updateTempPlot(heatsimFigure,model.HS.T);
          h_title.String = ['Temperature evolution, t = ' num2str(updatesTimeVector(updateIdx+1),'%#.2g') ' s'];
          drawnow;

          if model.HS.makeMovie
            movieFrames(movieFrameidx) = getframe(heatsimFigure);
            movieFrameidx = movieFrameidx + 1;
          end
        end

        if rem(updateIdx,model.HS.mediaPropRecalcPeriod) == 0 % If this update is one where we have to check for recalculating the optical or thermal properties
          [model,heatSimParameters] = heatSimParamsRecalc(model,heatSimParameters); % Recalculate MC and FMC, if necessary, and get new thermal properties and a new media matrix
        end
      end
      if ~model.HS.silentMode; fprintf(1,repmat('\b',1,nCharsOnLine)); end
    end
  end

  %% End relaxation phase
  if ~model.HS.silentMode && model.HS.durationEnd ~= 0
    updateIdx = model.HS.nPulses*(nUpdatesOn+nUpdatesOff); % Index of the update
    if updateIdx > 0
      timeRemainingEstimate = toc/updateIdx*((nUpdatesOn+nUpdatesOff)*model.HS.nPulses+nUpdatesEnd - updateIdx);
      nCharsOnLine = fprintf(1,'[%02.0f:%02.0f:%02.0f] Diffusing heat in end relaxation phase... 0/%d',floor(timeRemainingEstimate/3600),floor(rem(timeRemainingEstimate/60,60)),floor(rem(timeRemainingEstimate,60)),nUpdatesEnd);
    else
      nCharsOnLine = fprintf(1,'[00:00:00] Diffusing heat in end relaxation phase... 0/%d',nUpdatesEnd);
    end
    drawnow;
  end

  heatSimParameters.lightsOn = false;
  heatSimParameters.steps = nTsPerUpdateEnd;
  heatSimParameters.dt = dtEnd;
  for i = 1:nUpdatesEnd
    if model.HS.useGPU
      [model.HS.T,model.HS.Omega,newSensorTemps,newMaxMediaTemps] = finiteElementHeatPropagator_CUDA(model.HS.T,model.HS.Omega,heatSimParameters);
    else
      if forceSingleThreaded % Multithreading doesn't work properly on Linux older than R2020b for some reason
        [model.HS.T,model.HS.Omega,newSensorTemps,newMaxMediaTemps] = finiteElementHeatPropagator_singlethreaded(model.HS.T,model.HS.Omega,heatSimParameters);
      else
        [model.HS.T,model.HS.Omega,newSensorTemps,newMaxMediaTemps] = finiteElementHeatPropagator(model.HS.T,model.HS.Omega,heatSimParameters);
      end
    end

    SMidx = 1; % sub-medium index
    for iM=1:nM
      model.HS.maxMediaTemps(iM) = max([model.HS.maxMediaTemps(iM), double(newMaxMediaTemps(SMidx:(SMidx + mP_fHtrim(iM).nBins - 1)))]);
      SMidx = SMidx + mP_fHtrim(iM).nBins;
    end

    if isempty(model.HS.sensorTemps)
      model.HS.sensorTemps = newSensorTemps;
    else
      model.HS.sensorTemps = [model.HS.sensorTemps(:,1:end-1) newSensorTemps];
    end

    updateIdx = i + model.HS.nPulses*(nUpdatesOn+nUpdatesOff); % Index of the update
    if ~model.HS.silentMode
      fprintf(repmat('\b',1,nCharsOnLine));
      timeRemainingEstimate = toc/updateIdx*((nUpdatesOn+nUpdatesOff)*model.HS.nPulses+nUpdatesEnd - updateIdx);
      nCharsOnLine = fprintf(1,'[%02.0f:%02.0f:%02.0f] Diffusing heat in end relaxation phase... %d/%d',floor(timeRemainingEstimate/3600),floor(rem(timeRemainingEstimate/60,60)),floor(rem(timeRemainingEstimate,60)),i,nUpdatesEnd);
      updateTempPlot(heatsimFigure,model.HS.T);
      h_title.String = ['Temperature evolution, t = ' num2str(updatesTimeVector(updateIdx+1),'%#.2g') ' s'];
      drawnow;

      if model.HS.makeMovie
        movieFrames(movieFrameidx) = getframe(heatsimFigure);
        movieFrameidx = movieFrameidx + 1;
      end
    end

    if rem(updateIdx,model.HS.mediaPropRecalcPeriod) == 0 % If this update is one where we have to check for recalculating the optical or thermal properties
      [model,heatSimParameters] = heatSimParamsRecalc(model,heatSimParameters); % Recalculate MC and FMC, if necessary, and get new thermal properties and a new media matrix
    end
  end
  if ~model.HS.silentMode
    if model.HS.durationEnd ~= 0
      fprintf(1,repmat('\b',1,nCharsOnLine));
    end
    toc;
    set(heatsimFigure,'WindowStyle','Docked');
  end

  clear heatSimParameters;

  for iM=1:nM
    if model.HS.maxMediaTemps(iM) > 100
      fprintf("\nWarning: %s reaches a temperature above 100°C. If this medium contains\nwater, it will likely undergo changes not modeled in this code.\n\n", mP_fHtrim(iM).name);
    end
  end

  %% Finalize and write movie
  if ~model.HS.silentMode && model.HS.makeMovie
    if isstruct(model.HS.movieFrames(1)) % If there exists prior heat simulation movieFrames, append to them
      movieFrames = [model.HS.movieFrames movieFrames];
    end

    if model.HS.deferMovieWrite
      model.HS.movieFrames = movieFrames;
    else
      [resy,resx,~] = size(movieFrames(1).cdata);
      for idx = 1:length(movieFrames)
        movieFrames(idx).cdata = movieFrames(idx).cdata(1:end-rem(resy,2),1:end-rem(resx,2),:); % Crop frames by one pixel if x or y size is odd
      end
      movieFrames = [repmat(movieFrames(1),1,30) movieFrames(1:end) repmat(movieFrames(end),1,30)];
      writerObj = VideoWriter([cd '/temp.avi'],'Uncompressed AVI');
      fprintf('Writing video. This may take a few minutes and require a lot of temporary disk space...\n');
      open(writerObj);
      writeVideo(writerObj,movieFrames);
      close(writerObj);
      if contains(computer, 'MAC') % macOS operating system
        [~,~] = system([cd '/+MCmatlab/x264/x264_macOS -o "' cd '/' model.name '_heatSimoutput.mkv" "' cd '/temp.avi"']);
      elseif contains(computer, 'WIN') % Windows operating system
        [~,~] = system([cd '\+MCmatlab\x264\x264_win64.exe -o "' cd '\' model.name '_heatSimoutput.mkv" "' cd '\temp.avi"']);
      else % Linux operating system
        [~,~] = system([cd '/+MCmatlab/x264/x264_linux -o "' cd '/' model.name '_heatSimoutput.mkv" "' cd '/temp.avi"']);
      end
      delete([cd '/temp.avi']);
      fprintf('\b Done\n');
    end
  end
end

function [dTdtperdeltaT , dTdt_abs] = calcdTdtArrays(model)
  mP = model.HS.mediaProperties;

  VHC = single(mP.VHC); % Volumetric heat capacity array. Element i is the VHC of medium i in the split mediaProperties list.
  TC = single(mP.TC); % Thermal conductivity array. Element i is the TC of medium i in the split mediaProperties list.

  nSM = numel(mP.VHC);

  %% Calculate proportionality between voxel-to-voxel temperature difference DeltaT and time step temperature change dT
  TC_eff = 2*(TC*TC.')./(ones(nSM)*diag(TC)+diag(TC)*ones(nSM)); % Same as TC_eff(i,j) = 2*TC_red(i)*TC_red(j)/(TC_red(i)+TC_red(j)) but without for loops
  TC_eff(isnan(TC_eff)) = 0; % Neighboring insulating voxels return NaN but should just be 0

  dTdtperdeltaT  = cat(3,TC_eff./(diag(VHC)*ones(nSM))/model.G.dx^2,...
    TC_eff./(diag(VHC)*ones(nSM))/model.G.dy^2,...
    TC_eff./(diag(VHC)*ones(nSM))/model.G.dz^2); % [deg C /(s*deg C)] Time derivative of voxel temperature per voxel-to-voxel temperature difference. Third dimension corresponds to the different directions of heat diffusion (x, y and z)

  %% Calculate temperature change due to absorbed heat per time step
  % dTdt_abs [deg C/s] is the time derivative of voxel temperature due to absorption
  if ~isscalar(model.FMC.NFR)
    dTdt_abs = (sum(model.MC.NA,4) ... % Include absorption of excitation light, summed over all wavelengths...
      - model.FMC.sourceDistribution ... % but don't count the energy locally re-emitted as fluorescence...
      + sum(model.FMC.NA,4))*model.MC.P./VHC(model.HS.M); % but do include the final absorption of fluorescence light
  else
    dTdt_abs = sum(model.MC.NA,4)*model.MC.P./VHC(model.HS.M);
  end
end

function mindtmax = calcdtmax(M,TC,VHC,dx,dy,dz)
  TC  = single(TC(M));
  effectiveTCx = zeros(size(M) + [1 0 0],'single');
  effectiveTCy = zeros(size(M) + [0 1 0],'single');
  effectiveTCz = zeros(size(M) + [0 0 1],'single');
  effectiveTCx(2:end-1,:,:) = 2*TC(1:end-1,:,:).*TC(2:end,:,:)./(TC(1:end-1,:,:)+TC(2:end,:,:));
  effectiveTCy(:,2:end-1,:) = 2*TC(:,1:end-1,:).*TC(:,2:end,:)./(TC(:,1:end-1,:)+TC(:,2:end,:));
  effectiveTCz(:,:,2:end-1) = 2*TC(:,:,1:end-1).*TC(:,:,2:end)./(TC(:,:,1:end-1)+TC(:,:,2:end));
  clear TC

  effectiveTCx(isnan(effectiveTCx)) = 0; % Neighboring insulating voxels would return NaN but should just be 0
  effectiveTCy(isnan(effectiveTCy)) = 0;
  effectiveTCz(isnan(effectiveTCz)) = 0;

  individual_dtmax = VHC(M)./(effectiveTCx(1:end-1,:,:)/dx^2 + effectiveTCx(2:end,:,:)/dx^2 ...
    + effectiveTCy(:,1:end-1,:)/dy^2 + effectiveTCy(:,2:end,:)/dy^2 ...
    + effectiveTCz(:,:,1:end-1)/dz^2 + effectiveTCz(:,:,2:end)/dz^2);
  mindtmax = double(min(individual_dtmax(:)));
end

function [model,heatSimParameters] = heatSimParamsRecalc(model,heatSimParameters)
  MCsilent = model.MC.silentMode;
  FMCsilent = model.FMC.silentMode;

  if any(model.G.optTdependent(:) | model.G.optFDdependent(:)) % If there is a heat simulation dependence in the MC or FMC step...
    model.MC.silentMode = true;
    fprintf('(Recalculating Monte Carlo)');
    model = runMonteCarlo(model); % then run the MC again
    fprintf(repmat('\b',1,27));
    model.MC.silentMode = MCsilent;
    if ~isnan(model.FMC.simulationTime)
      model.FMC.silentMode = true;
      fprintf('(Recalculating Fluorescence Monte Carlo)');
      model = runMonteCarlo(model,'fluorescence'); % and if we have to simulate fluorescence, it needs to be recalculated too
      fprintf(repmat('\b',1,40));
      model.FMC.silentMode = FMCsilent;
    end
    model = getThermalMediaProperties(model); % Then get the new thermal media properties, including a new split media matrix
  elseif any(model.G.thmTdependent(:) | model.G.thmFDdependent(:)) % Otherwise if only the heat sim step (the thermal media properties) is temperature or FD dependent, no need to re-run any MC but the thermal media properties must be updated, including a new split media matrix.
    model = getThermalMediaProperties(model);
  end
  heatSimParameters.M = model.HS.M-1; % The new media matrix must be loaded into heatSimParameters
  [dTdtperdeltaT , dTdt_abs] = calcdTdtArrays(model);
  heatSimParameters.dTdtperdeltaT = dTdtperdeltaT; % Same for dTdtperdeltaT
  heatSimParameters.dTdt_abs = dTdt_abs; % and dTdt_abs
end

function [li , w] = convertTempSensorPos(tempSensorPositions,x,y,z)
  % To interpolate the temperature sensor temperatures during heat
  % simulation, the C script is going to need temperature values from the 8
  % surrounding voxels, calculated in the MEX file with this formula:
  % T_sensor =
  %     (1-wx)*(1-wy)*(1-wz)*T(ix  ,iy  ,iz  ) + (1-wx)*(1-wy)*wz*T(ix  ,iy  ,iz+1) +
  %     (1-wx)*   wy *(1-wz)*T(ix  ,iy+1,iz  ) + (1-wx)*   wy *wz*T(ix  ,iy+1,iz+1) +
  %        wx *(1-wy)*(1-wz)*T(ix+1,iy  ,iz  ) +    wx *(1-wy)*wz*T(ix+1,iy  ,iz+1) +
  %        wx *   wy *(1-wz)*T(ix+1,iy+1,iz  ) +    wx *   wy *wz*T(ix+1,iy+1,iz+1);
  %
  % This function calculates the linear index corresponding to the subscript
  % indices (ix,iy,iz) and the weights wx, wy and wz. The index is calculated
  % in C-notation, that is, starting from 0.

  nx = length(x);
  ny = length(y);
  nz = length(z);
  dx = x(2)-x(1);
  dy = y(2)-y(1);
  dz = z(2)-z(1);

  nS = size(tempSensorPositions,1); % Number of sensors
  li = zeros(nS,1); % Linear index of interpolation group corner voxel
  w = zeros(nS,3); % Weights in the x, y and z directions used in the interpolation

  for i=1:nS
    relposx = tempSensorPositions(i,1) - x(1);
    if relposx < 0
      ix = 0;
      w(i,1) = 0;
    elseif relposx < x(end) - x(1)
      ix = floor(relposx/dx);
      w(i,1) = mod(relposx/dx,1);
    else
      ix = nx-2;
      w(i,1) = 1;
    end

    relposy = tempSensorPositions(i,2) - y(1);
    if relposy < 0
      iy = 0;
      w(i,2) = 0;
    elseif relposy < y(end) - y(1)
      iy = floor(relposy/dy);
      w(i,2) = mod(relposy/dy,1);
    else
      iy = ny-2;
      w(i,2) = 1;
    end

    relposz = tempSensorPositions(i,3) - z(1);
    if relposz < 0
      iz = 0;
      w(i,3) = 0;
    elseif relposz < z(end) - z(1)
      iz = floor(relposz/dz);
      w(i,3) = mod(relposz/dz,1);
    else
      iz = nz-2;
      w(i,3) = 1;
    end

    li(i) = ix    + ...
      iy*nx + ...
      iz*nx*ny;
  end
end

function updateTempPlot(h_f,M)
  h_f.UserData(1:end-1,1:end-1,1:end-1) = single(M);
  h_f.UserData(end,:,:) = h_f.UserData(end-1,:,:);
  h_f.UserData(:,end,:) = h_f.UserData(:,end-1,:);
  h_f.UserData(:,:,end) = h_f.UserData(:,:,end-1);

  iC = 1; % Index of child
  while iC<length(h_f.Children) % Find out which graphics element is the handle for the checkbox
    if isprop(h_f.Children(iC),'Style') && strcmp(h_f.Children(iC).Style,'checkbox')
      v = h_f.Children(iC).Callback{2};
      break;
    else
      iC = iC+1;
    end
  end
  xsi = floor(get(v.h_sliders(1),'Value'));
  ysi = floor(get(v.h_sliders(2),'Value'));
  zsi = floor(get(v.h_sliders(3),'Value'));
  set(v.h_surfxslice,'CData',squeeze(h_f.UserData(xsi,:,:)));
  set(v.h_surfyslice,'CData',squeeze(h_f.UserData(:,ysi,:)));
  set(v.h_surfzslice,'CData',squeeze(h_f.UserData(:,:,zsi)));
  set(v.h_surfxback,'CData',squeeze(h_f.UserData(v.xbi,:,:)));
  set(v.h_surfyback,'CData',squeeze(h_f.UserData(:,v.ybi,:)));
  set(v.h_surfzback,'CData',squeeze(h_f.UserData(:,:,v.zbi)));
  drawnow;
end