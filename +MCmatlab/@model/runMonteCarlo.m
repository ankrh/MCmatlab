function model = runMonteCarlo(model,varargin)
%   Requires
%       MCmatlab.mex (architecture specific)
%
%   See also defineGeometry, plotMCmatlab, runMonteCarloFluorescence, simulateHeatDistribution

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

G = model.G;

if any(strcmp(varargin,'fluorescence'))
  simType = 2;
else
  simType = 1;
end

checkMCinputFields(model,simType);
model = getMediaProperties_funcHandles(model,simType); % Calls the mediaPropertiesFunc and converts those fields that are specified as char array formulas into function handles taking intensity and temperature as input

%% Check for GPU compatibility if needed
if (simType == 1 && model.MC.useGPU) || (simType == 2 && model.FMC.useGPU)
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

%% Choose what to use for initial temperature, fractional damage and fluence rate
if ((simType == 2 && model.FMC.Tdependent) || (simType == 1 && model.MC.Tdependent))
  if ~isnan(model.HS.T(1))
    T = model.HS.T;
  elseif isnan(model.HS.Tinitial(1))
    error('Error: Parameters are temperature dependent but no initial temperature is provided');
  else
    if numel(model.HS.Tinitial) == 1
      T = model.HS.Tinitial*ones(size(G.M_raw),'single');
    else
      T = single(model.HS.Tinitial);
    end
  end
else
  T = zeros(size(G.M_raw),'single');
end

if ~isnan(model.HS.Omega(1))
  FD = 1 - exp(-model.HS.Omega); % Fractional damage of molecules/cells
else
  FD = zeros(size(G.M_raw),'single');
end

if (simType == 2 && model.FMC.FRdependent) || (simType == 1 && model.MC.FRdependent)
  if isnan(model.MC.P)
    error('Error: Media properties are fluence rate dependent, but no power (model.MC.P) was specified');
  elseif ~isnan(model.MC.NFR(1))
    model.MC.FR = model.MC.P*model.MC.NFR; % If the MC simulations have been run before, use the most recent NFR
  elseif ~isnan(model.MC.FRinitial(1))
    model.MC.FR = model.MC.FRinitial; % Otherwise if it exists, use FRinitial as the initial guess
  else
    model.MC.FR = zeros(size(G.M_raw)); % Otherwise just start with no light
  end
end

%% Call Monte Carlo C code (MEX file) to get fluence rate distribution, possibly inside an iterative loop in the case of fluence rate dependent excitation light
iterate = simType == 1 && model.MC.FRdependent;

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
    if (simType == 1 && model.MC.useGPU) || (simType == 2 && model.FMC.useGPU)
      model = MCmatlab_CUDA(model,simType);
    else
      model = MCmatlab(model,simType);
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
    mua_vec = [model.MC.mediaProperties.mua];
    PY_3d = NaN(size(G.M_raw)); % 3D distribution of power yields
    for i=1:length(model.MC.mediaProperties_funcHandles)
      if isa(model.MC.mediaProperties_funcHandles(i).PY,'function_handle')
        if isnan(model.MC.FR(1))
          PY_3d(G.M_raw == i) = model.MC.mediaProperties_funcHandles(i).PY(zeros(size(G.M_raw)),T(G.M_raw == i),FD(G.M_raw(:) == i));
        else
          PY_3d(G.M_raw == i) = model.MC.mediaProperties_funcHandles(i).PY(model.MC.FR(G.M_raw == i),T(G.M_raw == i),FD(G.M_raw(:) == i));
        end
      elseif isa(model.MC.mediaProperties_funcHandles(i).QY,'function_handle')
        if isnan(model.MC.FR(1))
          PY_3d(G.M_raw == i) = model.MC.mediaProperties_funcHandles(i).QY(zeros(size(G.M_raw)),T(G.M_raw == i),FD(G.M_raw(:) == i))*model.MC.wavelength/model.FMC.wavelength;
        else
          PY_3d(G.M_raw == i) = model.MC.mediaProperties_funcHandles(i).QY(model.MC.FR(G.M_raw == i),T(G.M_raw == i),FD(G.M_raw(:) == i))*model.MC.wavelength/model.FMC.wavelength;
        end
      elseif ~isnan(model.MC.mediaProperties_funcHandles(i).PY)
        PY_3d(G.M_raw == i) = model.MC.mediaProperties_funcHandles(i).PY;
      elseif ~isnan(model.MC.mediaProperties_funcHandles(i).QY)
        PY_3d(G.M_raw == i) = model.MC.mediaProperties_funcHandles(i).QY*model.MC.wavelength/model.FMC.wavelength;
      end
    end
    PY_3d(isnan(PY_3d)) = 0; % Set any media that have no PY specified to be non-fluorescing
    if any(isinf(PY_3d(:)) | PY_3d(:) < 0)
      error('Error: Invalid fluorescence power or quantum yields');
    end
    model.FMC.beam.sourceDistribution = PY_3d.*mua_vec(model.MC.M).*model.MC.NFR.^model.FMC.fluorescenceOrder; % [W/cm^3/W.incident]
    clear PY_3d
    if max(model.FMC.beam.sourceDistribution(:)) == 0; error('Error: No fluorescence emitters'); end
  end
  if (simType == 1 && model.MC.useGPU) || (simType == 2 && model.FMC.useGPU)
    model = MCmatlab_CUDA(model,simType);
  else
    model = MCmatlab(model,simType);
  end
end

% Add positions of the centers of the pixels in the light collector image
if simType == 2
  if model.FMC.useLightCollector && model.FMC.LC.res > 1
    model.FMC.LC.X = linspace(model.FMC.LC.fieldSize*(1/model.FMC.LC.res-1),model.FMC.LC.fieldSize*(1-1/model.FMC.LC.res),model.FMC.LC.res)/2;
    model.FMC.LC.Y = model.FMC.LC.X;
  end

  % Add angles of the centers of the far field pixels
  if model.FMC.farFieldRes
    model.FMC.farFieldTheta = linspace(pi/model.FMC.farFieldRes/2,pi-pi/model.FMC.farFieldRes/2,model.FMC.farFieldRes);
    model.FMC.farFieldPhi   = linspace(-pi+(2*pi)/model.FMC.farFieldRes/2,pi-(2*pi)/model.FMC.farFieldRes/2,model.FMC.farFieldRes);
  end
else
  if model.MC.useLightCollector && model.MC.LC.res > 1
    model.MC.LC.X = linspace(model.MC.LC.fieldSize*(1/model.MC.LC.res-1),model.MC.LC.fieldSize*(1-1/model.MC.LC.res),model.MC.LC.res)/2;
    model.MC.LC.Y = model.MC.LC.X;
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

if isnan(MCorFMC.wavelength)
  error('Error: No wavelength defined');
end
if ~MCorFMC.calcNFR && ~MCorFMC.useLightCollector
  error('Error: calcNFR is false, but no light collector is defined');
end
if MCorFMC.calcNFRdet && ~MCorFMC.useLightCollector
  error('Error: calcNFRdet is true, but no light collector is defined');
end
if MCorFMC.farFieldRes && MCorFMC.boundaryType == 0
  error('Error: If boundaryType == 0, no photons can escape to be registered in the far field. Set farFieldRes to zero or change boundaryType.');
end

if simType == 2
  if isnan(model.MC.NFR(1))
    error('Error: NFR matrix not calculated for excitation light');
  end
else
  
  if isnan(MCorFMC.beam.beamType)
    error('Error: No beamType defined');
  end
  if (MCorFMC.beam.beamType == 3 || MCorFMC.beam.beamType == 4) && (isnan(MCorFMC.beam.NF.radialWidth) || isnan(MCorFMC.beam.FF.radialWidth))
    error('Error: beam.NF.radialWidth and beam.FF.radialWidth must both be specified when beamType is 3 or 4');
  end
  if MCorFMC.beam.beamType == 4 && (isnan(MCorFMC.beam.NF.radialDistr(1)) || isnan(MCorFMC.beam.FF.radialDistr(1)))
    error('Error: beam.NF.radialDistr and beam.FF.radialDistr must both be specified when beamType is 4');
  end
  if MCorFMC.beam.beamType == 5
    if isnan(MCorFMC.beam.NF.XDistr(1))
      error('Error: beam.NF.XDistr must be specified when beamType is 5');
    end
    if isnan(MCorFMC.beam.NF.XWidth)
      error('Error: beam.NF.Xwidth must be specified when beamType is 5');
    end
    if isnan(MCorFMC.beam.NF.YDistr(1))
      error('Error: beam.NF.YDistr must be specified when beamType is 5');
    end
    if isnan(MCorFMC.beam.NF.YWidth)
      error('Error: beam.NF.YWidth must be specified when beamType is 5');
    end
    if isnan(MCorFMC.beam.FF.XDistr(1))
      error('Error: beam.FF.XDistr must be specified when beamType is 5');
    end
    if isnan(MCorFMC.beam.FF.XWidth)
      error('Error: beam.FF.Xwidth must be specified when beamType is 5');
    end
    if isnan(MCorFMC.beam.FF.YDistr(1))
      error('Error: beam.FF.YDistr must be specified when beamType is 5');
    end
    if isnan(MCorFMC.beam.FF.YWidth)
      error('Error: beam.FF.YWidth must be specified when beamType is 5');
    end
    if xor(MCorFMC.beam.FF.XDistr == 2,MCorFMC.beam.FF.YDistr == 2)
      error('Error: beam.FF.XDistr and beam.FF.YDistr must either both be set to cosine (Lambertian), or neither');
    end
  end
  if size(MCorFMC.beam.NF.radialDistr,1) > 1 || size(MCorFMC.beam.NF.XDistr,1) > 1 || size(MCorFMC.beam.NF.YDistr,1) > 1 ||...
     size(MCorFMC.beam.FF.radialDistr,1) > 1 || size(MCorFMC.beam.FF.XDistr,1) > 1 || size(MCorFMC.beam.FF.YDistr,1) > 1
    error('Error: Beam distribution functions must be scalars or row-vectors, not column-vectors');
  end
  if (numel(MCorFMC.beam.NF.radialDistr) > 1 && numel(MCorFMC.beam.NF.radialDistr) < 1000) ||...
     (numel(MCorFMC.beam.NF.XDistr)      > 1 && numel(MCorFMC.beam.NF.XDistr)      < 1000) ||...
     (numel(MCorFMC.beam.NF.YDistr)      > 1 && numel(MCorFMC.beam.NF.YDistr)      < 1000) ||...
     (numel(MCorFMC.beam.FF.radialDistr) > 1 && numel(MCorFMC.beam.FF.radialDistr) < 1000) ||...
     (numel(MCorFMC.beam.FF.XDistr)      > 1 && numel(MCorFMC.beam.FF.XDistr)      < 1000) ||...
     (numel(MCorFMC.beam.FF.YDistr)      > 1 && numel(MCorFMC.beam.FF.YDistr)      < 1000)
    error('Error: Beam definition distributions must have 1000 elements (or more, if necessary)');
  end
  if any([MCorFMC.beam.NF.radialDistr , MCorFMC.beam.NF.XDistr , MCorFMC.beam.NF.YDistr , ...
          MCorFMC.beam.FF.radialDistr , MCorFMC.beam.FF.XDistr , MCorFMC.beam.FF.YDistr] < 0)
    error('Error: Beam definition distributions may not contain negative numbers');
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
      error('Error: LC.res must be 1 when LC.f is Inf');
    end
  end

  if (abs(xLCC)               < G.nx*G.dx/2 && ...
      abs(yLCC)               < G.ny*G.dy/2 && ...
      abs(zLCC - G.nz*G.dz/2) < G.nz*G.dz/2)
    error('Error: Light collector center (%.4f,%.4f,%.4f) is inside cuboid',xLCC,yLCC,zLCC);
  end
end
end