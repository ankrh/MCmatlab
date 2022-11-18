classdef geometry
  %GEOMETRY This class contains all properties related to the geometry of
  %an MCmatlab.model.
  %   This class defines the properties of a geometry to be used in a
  %   monteCarloSimulation, fluorescenceMonteCarloSimulation or
  %   heatSimulation.

  properties
    silentMode (1,1) logical = false % Disables command window text and progress indication
    nx (1,1) double {mustBeInteger, mustBeGreaterThan(nx,1)} = 100 % Number of bins in the x direction
    ny (1,1) double {mustBeInteger, mustBeGreaterThan(ny,1)} = 100 % Number of bins in the y direction
    nz (1,1) double {mustBeInteger, mustBeGreaterThan(nz,1)} = 100 % Number of bins in the z direction
    Lx (1,1) double {mustBeFinite, mustBePositive} = 1 % [cm] x size of simulation cuboid
    Ly (1,1) double {mustBeFinite, mustBePositive} = 1 % [cm] y size of simulation cuboid
    Lz (1,1) double {mustBeFinite, mustBePositive} = 1 % [cm] z size of simulation cuboid
    mediaPropertiesFunc (1,1) function_handle {mustHave1Input1OutputArgs} = @defaultMediaPropertiesFunc % Media properties defined as a function at the end of the model file
    mediaPropParams cell = {} % Cell array containing any additional parameters to be passed to the getMediaProperties function
    geomFunc (1,1) function_handle {mustHave4Input1OutputArgs} = @defaultGeomFunc % Function to use for defining the distribution of media in the cuboid. Defined at the end of the model file.
    geomFuncParams cell = {} % Cell array containing any additional parameters to pass into the geometry function, such as media depths, inhomogeneity positions, radii etc.

    %% Calculated
    M_raw uint8
    FRdependent (:,5) logical
    optTdependent (:,5) logical
    optFDdependent (:,5) logical
    thmTdependent (:,2) logical
    thmFDdependent (:,2) logical
  end

  properties (Hidden)
    needsRecalculation (1,1) logical = true
  end

  properties (Dependent)
    dx
    dy
    dz
    x
    y
    z
  end

  methods
    function obj = updateGeometry(obj)
      if obj.needsRecalculation
        [X,Y,Z] = ndgrid(single(obj.x),single(obj.y),single(obj.z)); % The single data type is used to conserve memory
        M_rawtemp = obj.geomFunc(X,Y,Z,obj.geomFuncParams);
        if ~isequal(size(M_rawtemp),[obj.nx, obj.ny, obj.nz])
          error('Error: The M returned by the geometry function must be nx by ny by nz in size. It''s a good idea to use, e.g., M = ones(size(X)); in your geometry definition to set the size correctly.');
        end
        if ~isreal(M_rawtemp)
          error('Error: The M returned by the geometry function must not contain complex numbers.');
        end
        if any(M_rawtemp(:) <= 0)
          error('Error: The M returned by the geometry function must only contain positive numbers.');
        end
        if any(rem(M_rawtemp(:),1))
          error('Error: The M returned by the geometry function must only contain integers.');
        end
        if any(M_rawtemp(:)>255)
          error('Error: The M returned by the geometry function must not contain numbers over 255.');
        end
        obj.M_raw = M_rawtemp;

        mP_fH = obj.mediaPropertiesFunc(obj.mediaPropParams); % Unsplit media
        mP_fH = mP_fH(unique(obj.M_raw));
        nM = numel(mP_fH);
        obj.FRdependent = false(nM,5);
        obj.optTdependent  = false(nM,5);
        obj.optFDdependent = false(nM,5);
        testWavelengths = linspace(100e-9,10000e-9,10);
        testFRs_vec = [0 logspace(-20,20,9)];
        testTs_vec = linspace(-100,300,10);
        testFDs_vec = linspace(0,1,10);
        [testFRs,testTs,testFD] = ndgrid(testFRs_vec,testTs_vec,testFDs_vec);
        for iM = 1:nM
          for iL = 1:numel(testWavelengths)
            muas = mP_fH(iM).mua(testWavelengths(iL),testFRs,testTs,testFD);
            muss = mP_fH(iM).mus(testWavelengths(iL),testFRs,testTs,testFD);
            gs   = mP_fH(iM).g  (testWavelengths(iL),testFRs,testTs,testFD);
            QYs  = mP_fH(iM).QY (testWavelengths(iL),testFRs,testTs,testFD);
            ESs  = mP_fH(iM).ES (testWavelengths(iL),testFRs,testTs,testFD);
            obj.FRdependent(iM,1) = obj.FRdependent(iM,1) || any(~isequaln(max(muas,[],1), min(muas,[],1)));
            obj.FRdependent(iM,2) = obj.FRdependent(iM,2) || any(~isequaln(max(muss,[],1), min(muss,[],1)));
            obj.FRdependent(iM,3) = obj.FRdependent(iM,3) || any(~isequaln(max(gs,[],1), min(gs,[],1)));
            obj.FRdependent(iM,4) = obj.FRdependent(iM,4) || any(~isequaln(max(QYs,[],1), min(QYs,[],1)));
            obj.FRdependent(iM,5) = obj.FRdependent(iM,5) || any(~isequaln(max(ESs,[],1), min(ESs,[],1)));
            obj.optTdependent(iM,1) = obj.optTdependent(iM,1) || any(~isequaln(max(muas,[],2), min(muas,[],2)));
            obj.optTdependent(iM,2) = obj.optTdependent(iM,2) || any(~isequaln(max(muss,[],2), min(muss,[],2)));
            obj.optTdependent(iM,3) = obj.optTdependent(iM,3) || any(~isequaln(max(gs,[],2), min(gs,[],2)));
            obj.optTdependent(iM,4) = obj.optTdependent(iM,4) || any(~isequaln(max(QYs,[],2), min(QYs,[],2)));
            obj.optTdependent(iM,5) = obj.optTdependent(iM,5) || any(~isequaln(max(ESs,[],2), min(ESs,[],2)));
            obj.optFDdependent(iM,1) = obj.optFDdependent(iM,1) || any(~isequaln(max(muas,[],3), min(muas,[],3)));
            obj.optFDdependent(iM,2) = obj.optFDdependent(iM,2) || any(~isequaln(max(muss,[],3), min(muss,[],3)));
            obj.optFDdependent(iM,3) = obj.optFDdependent(iM,3) || any(~isequaln(max(gs,[],3), min(gs,[],3)));
            obj.optFDdependent(iM,4) = obj.optFDdependent(iM,4) || any(~isequaln(max(QYs,[],3), min(QYs,[],3)));
            obj.optFDdependent(iM,5) = obj.optFDdependent(iM,5) || any(~isequaln(max(ESs,[],3), min(ESs,[],3)));
          end
          VHCs = mP_fH(iM).VHC(testTs,testFD);
          TCs  = mP_fH(iM).TC (testTs,testFD);
          obj.thmTdependent(iM,1) = any(~isequaln(max(VHCs,[],1), min(VHCs,[],1)));
          obj.thmTdependent(iM,2) = any(~isequaln(max(TCs,[],1), min(TCs,[],1)));
          obj.thmFDdependent(iM,1) = any(~isequaln(max(VHCs,[],2), min(VHCs,[],2)));
          obj.thmFDdependent(iM,2) = any(~isequaln(max(TCs,[],2), min(TCs,[],2)));
        end
        obj.needsRecalculation = false;
      end
    end
    
    function obj = set.nx(obj,val); obj.nx = val; obj.needsRecalculation = true; end %#ok<*MCSUP> 
    function obj = set.ny(obj,val); obj.ny = val; obj.needsRecalculation = true; end
    function obj = set.nz(obj,val); obj.nz = val; obj.needsRecalculation = true; end
    function obj = set.Lx(obj,val); obj.Lx = val; obj.needsRecalculation = true; end
    function obj = set.Ly(obj,val); obj.Ly = val; obj.needsRecalculation = true; end
    function obj = set.Lz(obj,val); obj.Lz = val; obj.needsRecalculation = true; end
    function obj = set.mediaPropertiesFunc(obj,val); obj.mediaPropertiesFunc = val; obj.needsRecalculation = true; end
    function obj = set.mediaPropParams(obj,val); obj.mediaPropParams = val; obj.needsRecalculation = true; end
    function obj = set.geomFunc(obj,val); obj.geomFunc = val; obj.needsRecalculation = true; end
    function obj = set.geomFuncParams(obj,val); obj.geomFuncParams = val; obj.needsRecalculation = true; end

    function value = get.dx(obj)
      value = obj.Lx/obj.nx; % [cm] size of x bins
    end
    function value = get.dy(obj)
      value = obj.Ly/obj.ny; % [cm] size of y bins
    end
    function value = get.dz(obj)
      value = obj.Lz/obj.nz; % [cm] size of z bins
    end
    function value = get.x(obj)
      value = ((0:obj.nx-1)-(obj.nx-1)/2)*obj.dx; % [cm] x position of centers of voxels
    end
    function value = get.y(obj)
      value = ((0:obj.ny-1)-(obj.ny-1)/2)*obj.dy; % [cm] y position of centers of voxels
    end
    function value = get.z(obj)
      value = ((0:obj.nz-1)+1/2)*obj.dz; % [cm] z position of centers of voxels
    end
  end
end

function mustHave4Input1OutputArgs(f)
if nargin(f) ~= 4
  error('Error: The function must take exactly 4 input arguments (X,Y,Z,parameters).');
end
if nargout(f) ~= 1
  error('Error: The function must return exactly 1 output argument (M, the media array)');
end
end

function mustHave1Input1OutputArgs(f)
if nargin(f) ~= 1
  error('Error: The syntax for defining media properties has changed slightly. The function must now take exactly 1 input argument (media properties parameters). Any wavelength dependence must be specified as in the examples.');
end
if nargout(f) ~= 1
  error('Error: The function must return exactly 1 output argument (an array of mediumProperties objects, the collection of media optical and (optionally) thermal properties). See the examples on how to define this array.');
end
end

function M = defaultGeomFunc(X,~,~,~)
M = ones(size(X));
end

function mediaProperties = defaultMediaPropertiesFunc(~)
mediaProperties = MCmatlab.mediumProperties;
end