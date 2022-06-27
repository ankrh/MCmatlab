classdef model
  % MODEL Top-level MCmatlab model class
  %
  %   MCmatlab.model collects all objects (as properties) and methods
  %   required to define a full MCmatlab model.

  properties
    G   (1,1) MCmatlab.geometry
    MC  (1,1) MCmatlab.monteCarloSimulation
    FMC (1,1) MCmatlab.fluorescenceMonteCarloSimulation
    HS  (1,1) MCmatlab.heatSimulation
  end

  methods
    function obj = model()
      % MODEL Construct an instance of this class
      %
      %   Create a new, clean instance of MCmatlab.model

      obj.G = MCmatlab.geometry; % We need to initialize the G object because it is a handle object. The other objects (MC, FMC, HS etc.) are value objects and do not need to be initialized in this constructor.

      %% The following code is inspired by the version check method of the unrelated but excellent export_fig() function made by Yair Altman
      % Check for newer version (only once a day)
      persistent lastCheckTime
      if isempty(lastCheckTime) || now - lastCheckTime > 1
        try
          currentVersion = str2double(strsplit(fileread('version'),'.'));
          url = 'https://raw.githubusercontent.com/ankrh/MCmatlab/Release/version';
          latestVersion = str2int32(strsplit(webread(url),'.'));
          if latestVersion(1) >  currentVersion(1) || ...
            (latestVersion(1) == currentVersion(1) && latestVersion(2) > currentVersion(2))
            msg = sprintf(['You are using version %d.%d of MCmatlab. ' ...
              'A newer version (%d.%d) is available for download from GitHub.'], ...
              currentVersion(1), currentVersion(2), latestVersion(1), latestVersion(2));
            msg = hyperlink('https://github.com/ankrh/MCmatlab', 'GitHub', msg);
            warning('MCmatlab:version',msg);
          end
        catch
          % ignore
        end
      end
      lastCheckTime = now;
    end

  function obj = reset(obj, type)
    % RESET(obj, type) Reset one of the sub-objects of an
    % MCmatlab.model to its default values.
    %
    %   'obj' is an object of class MCmatlab.model.
    %   'type' is one of 'G', 'MC', 'FMC', or 'HS'.
    %
    %   Reset the properties of either geometry ('G'),
    %   monteCarloSimulation ('MC'), fluorescenceMonteCarloSimulation ('FMC'), or
    %   heatSimulation ('HS') in a MCmatlab.model to their defualt values.

    if nargin == 1
      obj = MCmatlab.model;
      return
    end

    switch type
      case "G"
        obj.G = MCmatlab.geometry;
      case "MC"
        obj.MC = MCmatlab.monteCarloSimulation;
      case "FMC"
        obj.FMC = MCmatlab.fluorescenceMonteCarloSimulation;
      case "HS"
        obj.HS = MCmatlab.heatSimulation;
      otherwise
        error("No valid type to reset specified (G, MC, FMC, HS)");
    end
  end


  function obj = clearMCmatlabModel(obj, type)
    warning('clearMCmatlabModel(obj, type) is deprecated, use reset(obj, type) instead.')

    reset(obj, type);
  end

  obj = runMonteCarlo(obj, varargin)

  function plot(obj, type)
    % PLOT(obj, type) Plot all calculated values
    %
    %   'obj' is an object of class MCmatlab.model.
    %   'type' is one of 'G', 'MC', 'FMC', or 'HS'.
    %
    %   Plots the calculated values of either geometry ('G'),
    %   monteCarloSimulation ('MC'), fluorescenceMonteCarloSimulation ('FMC'), or
    %   heatSimulation ('HS') in a MCmatlab.model.

    if nargin == 1
      % plot everything that is calculated
      if ~isnan(obj.G.nx); plotMCmatlabGeom(obj); end
      if ~isnan(obj.MC.nPhotons); plotMCmatlabMC(obj); end
      if ~isnan(obj.FMC.nPhotons); plotMCmatlabMC(obj, "fluorescence"); end
      if ~isnan(obj.HS.T); plotMCmatlabHeat(obj); end
      return
    end

    switch type
      case "FMC"
        plotMCmatlabMC(obj, "fluorescence")
      case "G"
        plotMCmatlabGeom(obj)
      case "HS"
        plotMCmatlabHeat(obj)
      case "MC"
        plotMCmatlabMC(obj)
      otherwise
        error("No valid plot type specified (G, MC, FMC, HS)");
    end

  end

  function plotMCmatlab(obj, varargin)
    warning("plotMCmatlab(obj, type) is deprecated, use plot(obj, type) instead.")

    if nargin == 1
      plotMCmatlabMC(obj)
      return
    end
    switch cell2mat(varargin)
      case "fluorescence"
        plotMCmatlabMC(obj, "fluorescence")
      case "geometry"
        plotMCmatlabGeom(obj)
      case "heatSimulation"
        plotMCmatlabHeat(obj)
    end
  end

  function obj = defineGeometry(obj)
    warning('Calling ''defineGeometry'' is no longer required. You can delete any calls to ''defineGeometry'' from your model files.')
  end

  obj = simulateHeatDistribution(obj)

  plotMCmatlabMC(obj, varargin)
  plotMCmatlabGeom(obj)
  plotMCmatlabHeat(obj)

end

methods (Static)
  function obj = loadobj(obj)
    warning('You have loaded an MCmatlab model from a mat-file. To enable working with the loaded model, you should now manually run the two lines from the original model file that set ''model.G.mediaPropertiesFunc'' and ''model.G.geomFunc''.')
  end
end

methods (Access = private)
  finiteElementHeatPropagator(T,Omega,heatSimParameters)
  finiteElementHeatPropagator_CUDA(T,Omega,heatSimParameters)

  getMediaProperties_funcHandles(model,simType)
  getOpticalMediaProperties(model,simType)
  getThermalMediaProperties(model)

  plotMediaProperties(nFig,model,simType)
  inferno(m)

  MCmatlab(model,simType)
  MCmatlab_CUDA(model,simType)
end

end

