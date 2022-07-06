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
    mediaPropertiesFunc (1,1) function_handle {mustHave2Input1OutputArgs} = @defaultMediaPropertiesFunc % Media properties defined as a function at the end of the model file
    mediaPropParams cell = {} % Cell array containing any additional parameters to be passed to the getMediaProperties function
    geomFunc (1,1) function_handle {mustHave4Input1OutputArgs} = @defaultGeomFunc % Function to use for defining the distribution of media in the cuboid. Defined at the end of the model file.
    geomFuncParams cell = {} % Cell array containing any additional parameters to pass into the geometry function, such as media depths, inhomogeneity positions, radii etc.
    M_raw uint8
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
    function obj = update_M_raw(obj)
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
    end
    
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

function mustHave2Input1OutputArgs(f)
if nargin(f) ~= 2
  error('Error: The function must take exactly 2 input arguments (wavelength,parameters).');
end
if nargout(f) ~= 1
  error('Error: The function must return exactly 1 output argument (mediaProperties, the struct of media optical and (optionally) thermal properties)');
end
end

function M = defaultGeomFunc(X,~,~,~)
M = ones(size(X));
end

function mediaProperties = defaultMediaPropertiesFunc(~,~)
mediaProperties = struct();
end