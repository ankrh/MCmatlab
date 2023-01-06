classdef model
  % MODEL Top-level MCmatlab model class
  %
  %   MCmatlab.model collects all objects (as properties) and methods
  %   required to define a full MCmatlab model.

  properties
    name (1,:) char
    G    (1,1) MCmatlab.geometry
    MC   (1,1) MCmatlab.monteCarloSimulation
    FMC  (1,1) MCmatlab.fluorescenceMonteCarloSimulation
    HS   (1,1) MCmatlab.heatSimulation
  end

  methods
    function obj = model()
      caller = dbstack(1);
      if isempty(caller) || contains(caller(1).name,'LiveEditorEvaluation')
        [~,obj.name,~] = fileparts(matlab.desktop.editor.getActiveFilename);
%         warning('Could not automatically determine the model name. Assuming the model name is that of the file currently open in the editor "%s". This might be caused by running the script with F9 or "Run Section". The model name is used only to determine which file to save output heat simulator videos to.',obj.name);
      else
        obj.name = caller(1).name;
      end

      addpath('./helperfuncs/');

      %% The following code is inspired by the version check method of the unrelated but excellent export_fig() function made by Yair Altman
      % Check for newer version (only once a day)
      persistent lastCheckTime
      if isempty(lastCheckTime) || now - lastCheckTime > 1
        try
          currentVersion = str2double(strsplit(fileread('version'),'.'));
          url = 'https://raw.githubusercontent.com/ankrh/MCmatlab/Release/version';
          latestVersion = str2double(strsplit(webread(url),'.'));
          newerVersion = false;
          for iNumber = 1:numel(currentVersion)
            if latestVersion(iNumber) > currentVersion(iNumber)
              newerVersion = true;
              break;
            elseif latestVersion(iNumber) < currentVersion(iNumber) % Something's wrong with the version number
              break;
            end
          end
          if newerVersion
            warning(['You are using version ' strrep(num2str(currentVersion),'  ','.') ' of MCmatlab. ' ...
              'A newer version (' strrep(num2str(latestVersion),'  ','.') ') is available for download from ' ...
              '<a href="matlab:web(''-browser'',''https://github.com/ankrh/MCmatlab'');">GitHub</a>.']);
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

    function obj = plot(obj, type)
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
        if ~isnan(obj.G.nx); obj = plotMCmatlabGeom(obj); end
        if ~isnan(obj.MC.nPhotons); obj = plotMCmatlabMC(obj); end
        if ~isnan(obj.FMC.nPhotons); obj = plotMCmatlabMC(obj, "fluorescence"); end
        if ~isnan(obj.HS.T); obj = plotMCmatlabHeat(obj); end
        return
      end

      if nargout == 0
        error(['Error: The syntax for calling plot has changed. The model object is now both an input and an output.' ...
          ' Use, e.g., "model = plot(model,''', type ,''')"']);
      end

      switch type
        case "FMC"
          obj = plotMCmatlabMC(obj, "fluorescence");
        case "G"
          obj = plotMCmatlabGeom(obj);
        case "HS"
          obj = plotMCmatlabHeat(obj);
        case "MC"
          obj = plotMCmatlabMC(obj);
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
    obj = plotMCmatlabMC(obj, varargin)
    obj = plotMCmatlabGeom(obj)
    obj = plotMCmatlabHeat(obj)
  end

  methods (Static)
    function obj = loadobj(obj)
      addpath('./helperfuncs/');
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

