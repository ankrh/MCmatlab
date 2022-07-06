function model = combineModels(modelArray,simType)
  arguments
    modelArray (1,:) MCmatlab.model
    simType (1,:) char
  end

  % Verify that G properties are identical
  Gprops = {'nx','ny','nz','Lx','Ly','Lz','mediaPropertiesFunc','mediaPropParams','geomFunc','geomFuncParams'};
  for iProp = 1:numel(Gprops)
    allEqual = true;
    for iModel = 2:numel(modelArray)
      allEqual = allEqual && isequaln(modelArray(iModel).G.(Gprops{iProp}),modelArray(1).G.(Gprops{iProp}));
    end
    if ~allEqual
      error('Error: The G %s properties are not identical, so the objects cannot be combined.',Gprops{iProp});
    end
  end

  % Verify that all MC sims have been run
  for iModel = 1:numel(modelArray)
    if isnan(modelArray(iModel).MC.nPhotons)
      error('Error: All input models must have had their MC simulation executed using runMonteCarlo() before combining.');
    end
  end

  % Verify that the relevant MC input properties are identical
  MCprops = {'calcNormalizedFluenceRate','calcNormalizedFluenceRate_detected','farFieldRes','matchedInterfaces','smoothingLengthScale','boundaryType','wavelength','useLightCollector','depositionCriteria'};
  for iProp = 1:numel(MCprops)
    allEqual = true;
    for iModel = 2:numel(modelArray)
      allEqual = allEqual && isequaln(modelArray(iModel).MC.(MCprops{iProp}),modelArray(1).MC.(MCprops{iProp}));
    end
    if ~allEqual
      error('Error: The MC %s properties are not identical, so the objects cannot be combined.',MCprops{iProp});
    end
  end

  % Check whether the light sources and powers are identical, in which case
  % we will weight the means by nPhotons. Otherwise we weight with the
  % powers
  lightSourcesAndPowersEqual = true;
  for iModel = 2:numel(modelArray)
    lightSourcesAndPowersEqual = lightSourcesAndPowersEqual && isequaln(modelArray(iModel).MC.LS,modelArray(1).MC.LS) ...
      && isequaln(modelArray(iModel).MC.P ,modelArray(1).MC.P );
  end

  % Verify that the relevant MC.LC input properties are identical
  LCprops = {'x','y','z','theta','phi','f','diam','fieldSize','NA','res','tStart','tEnd','nTimeBins'};
  for iProp = 1:numel(LCprops)
    allEqual = true;
    for iModel = 2:numel(modelArray)
      allEqual = allEqual && isequaln(modelArray(iModel).MC.LC.(LCprops{iProp}),modelArray(1).MC.LC.(LCprops{iProp}));
    end
    if ~allEqual
      error('Error: The light collector %s properties are not identical, so the objects cannot be combined.',LCprops{iProp});
    end
  end

  if strcmp(simType,'MC')
    %% Combine MC
    model = modelArray(1);
    reset(model,'FMC');
    reset(model,'HS');

    if ~lightSourcesAndPowersEqual
      if isnan(modelArray(1).MC.P)
        error('Error: If the light source definitions are not identical, all input models must have the MC.P (power) property specified. These do not need to be the correct absolute power values, but must have the correct relative power values.');
      end
      w1 = modelArray(1).MC.P; % Weight 1
    else
      w1 = modelArray(1).MC.nPhotons; % Weight 1
    end

    for iModel = 2:numel(modelArray)
      if ~lightSourcesAndPowersEqual
        if isnan(modelArray(iModel).MC.P)
          error('Error: If the light source definitions are not identical, all input models must have the MC.P (power) property specified. These do not need to be the correct absolute power values, but must have the correct relative power values.');
        end
        w = modelArray(iModel).MC.P; % Weight
      else
        w = modelArray(iModel).MC.nPhotons; % Weight
      end

      model.MC.simulationTime = model.MC.simulationTime + modelArray(iModel).MC.simulationTime;
      model.MC.nPhotons = model.MC.nPhotons + modelArray(iModel).MC.nPhotons;
      model.MC.examplePaths = cat(2,model.MC.examplePaths,modelArray(iModel).MC.examplePaths);
      if ~lightSourcesAndPowersEqual
        model.MC.P = model.MC.P + modelArray(iModel).MC.P;
      end

      model.MC.LC.image   = model.MC.LC.image   + w/w1*modelArray(iModel).MC.LC.image;
      model.MC.NFR        = model.MC.NFR        + w/w1*modelArray(iModel).MC.NFR;
      model.MC.farField   = model.MC.farField   + w/w1*modelArray(iModel).MC.farField;
      model.MC.NI_xpos    = model.MC.NI_xpos    + w/w1*modelArray(iModel).MC.NI_xpos;
      model.MC.NI_xneg    = model.MC.NI_xneg    + w/w1*modelArray(iModel).MC.NI_xneg;
      model.MC.NI_ypos    = model.MC.NI_ypos    + w/w1*modelArray(iModel).MC.NI_ypos;
      model.MC.NI_yneg    = model.MC.NI_yneg    + w/w1*modelArray(iModel).MC.NI_yneg;
      model.MC.NI_zpos    = model.MC.NI_zpos    + w/w1*modelArray(iModel).MC.NI_zpos;
      model.MC.NI_zneg    = model.MC.NI_zneg    + w/w1*modelArray(iModel).MC.NI_zneg;
    end

    if ~lightSourcesAndPowersEqual
      wTot = model.MC.P; % Weight total
      model.MC.LS = MCmatlab.lightSource.empty; % Since we now have a combination of light sources, we represent this by leaving the property empty
    else
      wTot = model.MC.nPhotons; % Weight total
    end

    % Renormalize outputs
    model.MC.LC.image   = model.MC.LC.image  *w1/wTot;
    model.MC.NFR        = model.MC.NFR       *w1/wTot;
    model.MC.farField   = model.MC.farField  *w1/wTot;
    model.MC.NI_xpos    = model.MC.NI_xpos   *w1/wTot;
    model.MC.NI_xneg    = model.MC.NI_xneg   *w1/wTot;
    model.MC.NI_ypos    = model.MC.NI_ypos   *w1/wTot;
    model.MC.NI_yneg    = model.MC.NI_yneg   *w1/wTot;
    model.MC.NI_zpos    = model.MC.NI_zpos   *w1/wTot;
    model.MC.NI_zneg    = model.MC.NI_zneg   *w1/wTot;
  elseif strcmp(simType,'FMC')
    %% Combine FMC
    % Verify that all MC objects are identical
    for iModel = 2:numel(modelArray)
      if ~isequaln(modelArray(iModel).MC,modelArray(1).MC)
        error('Error: The MC objects must all be identical for the FMC objects to be combined.');
      end
    end

    % Verify that all FMC sims have been run
    for iModel = 1:numel(modelArray)
      if isnan(modelArray(iModel).FMC.nPhotons)
        error('Error: All input models must have had their FMC simulation executed using runMonteCarlo() before combining.');
      end
    end

    % Verify that the relevant FMC input properties are identical
    FMCprops = {'calcNormalizedFluenceRate','calcNormalizedFluenceRate_detected','farFieldRes','matchedInterfaces','smoothingLengthScale','boundaryType','wavelength','useLightCollector'};
    for iProp = 1:numel(FMCprops)
      allEqual = true;
      for iModel = 2:numel(modelArray)
        allEqual = allEqual && isequaln(modelArray(iModel).FMC.(FMCprops{iProp}),modelArray(1).FMC.(FMCprops{iProp}));
      end
      if ~allEqual
        error('Error: The FMC %s properties are not identical, so the objects cannot be combined.',FMCprops{iProp});
      end
    end

    % Verify that the relevant FMC.LC input properties are identical
    LCprops = {'x','y','z','theta','phi','f','diam','fieldSize','NA','res','tStart','tEnd','nTimeBins'};
    for iProp = 1:numel(LCprops)
      allEqual = true;
      for iModel = 2:numel(modelArray)
        allEqual = allEqual && isequaln(modelArray(iModel).FMC.LC.(LCprops{iProp}),modelArray(1).FMC.LC.(LCprops{iProp}));
      end
      if ~allEqual
        error('Error: The FMC light collector %s properties are not identical, so the objects cannot be combined.',LCprops{iProp});
      end
    end

    model = modelArray(1);
    reset(model,'HS');

    w1 = modelArray(1).FMC.nPhotons; % Weight 1
    for iModel = 2:numel(modelArray)
      w = modelArray(iModel).FMC.nPhotons; % Weight

      model.FMC.simulationTime = model.FMC.simulationTime + modelArray(iModel).FMC.simulationTime;
      model.FMC.nPhotons = model.FMC.nPhotons + modelArray(iModel).FMC.nPhotons;
      model.FMC.examplePaths = cat(2,model.FMC.examplePaths,modelArray(iModel).FMC.examplePaths);

      model.FMC.LC.image   = model.FMC.LC.image   + w/w1*modelArray(iModel).FMC.LC.image;
      model.FMC.NFR        = model.FMC.NFR        + w/w1*modelArray(iModel).FMC.NFR;
      model.FMC.farField   = model.FMC.farField   + w/w1*modelArray(iModel).FMC.farField;
      model.FMC.NI_xpos    = model.FMC.NI_xpos    + w/w1*modelArray(iModel).FMC.NI_xpos;
      model.FMC.NI_xneg    = model.FMC.NI_xneg    + w/w1*modelArray(iModel).FMC.NI_xneg;
      model.FMC.NI_ypos    = model.FMC.NI_ypos    + w/w1*modelArray(iModel).FMC.NI_ypos;
      model.FMC.NI_yneg    = model.FMC.NI_yneg    + w/w1*modelArray(iModel).FMC.NI_yneg;
      model.FMC.NI_zpos    = model.FMC.NI_zpos    + w/w1*modelArray(iModel).FMC.NI_zpos;
      model.FMC.NI_zneg    = model.FMC.NI_zneg    + w/w1*modelArray(iModel).FMC.NI_zneg;
    end

    wTot = model.FMC.nPhotons; % Weight total

    % Renormalize outputs
    model.FMC.LC.image   = model.FMC.LC.image  *w1/wTot;
    model.FMC.NFR        = model.FMC.NFR       *w1/wTot;
    model.FMC.farField   = model.FMC.farField  *w1/wTot;
    model.FMC.NI_xpos    = model.FMC.NI_xpos   *w1/wTot;
    model.FMC.NI_xneg    = model.FMC.NI_xneg   *w1/wTot;
    model.FMC.NI_ypos    = model.FMC.NI_ypos   *w1/wTot;
    model.FMC.NI_yneg    = model.FMC.NI_yneg   *w1/wTot;
    model.FMC.NI_zpos    = model.FMC.NI_zpos   *w1/wTot;
    model.FMC.NI_zneg    = model.FMC.NI_zneg   *w1/wTot;
  end
end
