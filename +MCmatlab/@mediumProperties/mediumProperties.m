classdef mediumProperties
  properties
    name                (1,:) char {mustBeNonempty} = 'Unlabeled medium'
    mua                 (1,1) = @defaultnan_41func
    mus                 (1,1) = @defaultnan_41func
    g                   (1,1) = @defaultnan_41func
    customPhaseFunc     (1,1) = @defaultnan_21func
    n                   (1,1) = @default1_11func
    QY                  (1,1) = @default0_41func
    ES                  (1,1) = @default1_41func
    VHC                 (1,1) = @defaultnan_21func
    TC                  (1,1) = @defaultnan_21func
    E                   (1,1) double {mustBeNonnegative} = 0
    A                   (1,1) double {mustBeNonnegative} = 0
    nBins               (1,1) double {mustBeInteger,mustBePositive} = 1
  end

  properties (Hidden)
    PY
  end

  methods
    function obj = set.PY(~,~) %#ok<STOUT> 
      error('Error: Power yield has been deprecated. Please use quantum yield (QY) instead.');
    end
    function obj = set.mua(obj,x)
      if isa(x,'function_handle')
        if abs(nargout(x)) ~= 1
          error('Error: The function must return exactly one output argument (mua).')
        end
        if nargin(x) ~= 4 && nargin(x) ~= 1
          error('Error: The function must take either one input argument (wavelength) or four input arguments (wavelength, fluence rate (FR), temperature (T), fractional damage (FD)).');
        end
        if nargin(x) == 1
          obj.mua = @(a,~,~,~)(x(a)); % Assume the input is the wavelength
        else
          obj.mua = x;
        end
      elseif isnumeric(x)
        obj.mua = @(~,~,~,~)x;
      else
        error('Error: Unexpected input type.');
      end
    end

    function obj = set.mus(obj,x)
      if isa(x,'function_handle')
        if abs(nargout(x)) ~= 1
          error('Error: The function must return exactly one output argument (mus).')
        end
        if nargin(x) ~= 4 && nargin(x) ~= 1
          error('Error: The function must take either one input argument (wavelength) or four input arguments (wavelength, fluence rate (FR), temperature (T), fractional damage (FD)).');
        end
        if nargin(x) == 1
          obj.mus = @(a,~,~,~)(x(a)); % Assume the input is the wavelength
        else
          obj.mus = x;
        end
      elseif isnumeric(x)
        obj.mus = @(~,~,~,~)x;
      else
        error('Error: Unexpected input type.');
      end
    end

    function obj = set.g(obj,x)
      if isa(x,'function_handle')
        if abs(nargout(x)) ~= 1
          error('Error: The function must return exactly one output argument (g).')
        end
        if nargin(x) ~= 4 && nargin(x) ~= 1
          error('Error: The function must take either one input argument (wavelength) or four input arguments (wavelength, fluence rate (FR), temperature (T), fractional damage (FD)).');
        end
        if nargin(x) == 1
          obj.g = @(a,~,~,~)(x(a)); % Assume the input is the wavelength
        else
          obj.g = x;
        end
      elseif isnumeric(x)
        obj.g = @(~,~,~,~)x;
      else
        error('Error: Unexpected input type.');
      end
    end

    function obj = set.customPhaseFunc(obj,x)
      if isa(x,'function_handle')
        if abs(nargout(x)) ~= 1
          error('Error: The function must return exactly one output argument (phaseFunction).')
        end
        if nargin(x) ~= 2
          error('Error: The function must take exactly two input arguments (wavelength, theta).');
        end
        obj.customPhaseFunc = x;
      else
        error('Error: Unexpected input type.');
      end
    end
    
    function obj = set.n(obj,x)
      if isa(x,'function_handle')
        if abs(nargout(x)) ~= 1
          error('Error: The function must return exactly one output argument (n).')
        end
        if nargin(x) ~= 1
          error('Error: The function must take exactly one input argument (wavelength).');
        end
        if nargin(x) == 1
          obj.n = @(a,~,~,~)(x(a)); % Assume the input is the wavelength
        else
          obj.n = x;
        end
      elseif isnumeric(x)
        obj.n = @(~,~,~,~)x;
      else
        error('Error: Unexpected input type.');
      end
    end

    function obj = set.QY(obj,x)
      if isa(x,'function_handle')
        if abs(nargout(x)) ~= 1
          error('Error: The function must return exactly one output argument (QY, the fluorescence quantum yield as a function of excitation wavelength).')
        end
        if nargin(x) ~= 4 && nargin(x) ~= 1
          error('Error: The function must take either one input argument (wavelength) or four input arguments (excitation wavelength, fluence rate (FR), temperature (T), fractional damage (FD)).');
        end
        if nargin(x) == 1
          obj.QY = @(a,~,~,~)(x(a)); % Assume the input is the wavelength
        else
          obj.QY = x;
        end
      elseif isnumeric(x)
        obj.QY = @(~,~,~,~)x;
      else
        error('Error: Unexpected input type.');
      end
    end

    function obj = set.ES(obj,x)
      if isa(x,'function_handle')
        if abs(nargout(x)) ~= 1
          error('Error: The function must return exactly one output argument (emissionSpectrum).')
        end
        if nargin(x) ~= 4 && nargin(x) ~= 1
          error('Error: The function must take either one input argument (wavelength) or four input arguments (wavelength, fluence rate (FR), temperature (T), fractional damage (FD)).');
        end
        if nargin(x) == 1
          obj.ES = @(a,~,~,~)(x(a)); % Assume the input is the wavelength
        else
          obj.ES = x;
        end
      elseif isnumeric(x)
        obj.ES = @(~,~,~,~)x;
      else
        error('Error: Unexpected input type.');
      end
    end

    function obj = set.VHC(obj,x)
      if isa(x,'function_handle')
        if abs(nargout(x)) ~= 1
          error('Error: The function must return exactly one output argument (volumetricHeatCapacity).')
        end
        if nargin(x) ~= 2
          error('Error: The function must take exactly two input arguments (temperature (T), fractional damage (FD)).');
        end
        obj.VHC = x;
      elseif isnumeric(x)
        obj.VHC = @(~,~)x;
      else
        error('Error: Unexpected input type.');
      end
    end

    function obj = set.TC(obj,x)
      if isa(x,'function_handle')
        if abs(nargout(x)) ~= 1
          error('Error: The function must return exactly one output argument (thermalConductivity).')
        end
        if nargin(x) ~= 2
          error('Error: The function must take exactly two input arguments (temperature (T), fractional damage (FD)).');
        end
        obj.TC = x;
      elseif isnumeric(x)
        obj.TC = @(~,~)x;
      else
        error('Error: Unexpected input type.');
      end
    end
  end
end

function x = defaultnan_41func(~,~,~,~)
  x = NaN;
end

function x = default0_41func(~,~,~,~)
  x = 0;
end

function x = default1_41func(~,~,~,~)
  x = 1;
end

function x = default1_11func(~)
  x = 1;
end

function x = defaultnan_21func(~,~)
  x = NaN;
end
