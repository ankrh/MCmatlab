classdef sourceIntensityDistribution
  %   sourceIntensityDistribution Summary of this class goes here
  %   Detailed explanation goes here

  properties
    radialDistr (1,:) double {mustBe012OrNaNOrNonnegativeArray} = NaN % Radial near or far field distribution - 0: Top-hat, 1: Gaussian, 2: Cosine (Lambertian), Array: Custom. Doesn't need to be normalized.
    radialWidth (1,1) double {mustBeNonnegativeOrNaN} = NaN % [cm] Radial near field 1/e^2 radius if top-hat or Gaussian or half-width of the full distribution if custom
                                                            % [rad] Radial far field 1/e^2 half-angle if top-hat or Gaussian or half-angle of the full distribution if custom. For a diffraction limited Gaussian beam, this should be set to model.MC.wavelength*1e-9/(pi*model.MC.beam.NF.radialWidth*1e-2))
    XDistr (1,:) double {mustBe012OrNaNOrNonnegativeArray} = NaN % X near or far field distribution - 0: Top-hat, 1: Gaussian, 2: Cosine (Lambertian), Array: Custom. Doesn't need to be normalized.
    XWidth (1,1) double {mustBeNonnegativeOrNaN} = NaN % [cm] X near field 1/e^2 radius if top-hat or Gaussian or half-width of the full distribution if custom
                                                       % [rad] X far field 1/e^2 half-angle if top-hat or Gaussian or half-angle of the full distribution if custom
    YDistr (1,:) double {mustBe012OrNaNOrNonnegativeArray} = NaN % Y near or far field distribution - 0: Top-hat, 1: Gaussian, 2: Cosine (Lambertian), Array: Custom. Doesn't need to be normalized.
    YWidth (1,1) double {mustBeNonnegativeOrNaN} = NaN % [cm] Y near field 1/e^2 radius if top-hat or Gaussian or half-width of the full distribution if custom
                                                       % [rad] Y far field 1/e^2 half-angle if top-hat or Gaussian or half-angle of the full distribution if custom
  end

  methods

  end
end

function mustBe012OrNaNOrNonnegativeArray(x)
x = x(:);
if isscalar(x)
  if (x == 0) || (x == 1) || (x == 2) || isnan(x)
    % Input valid
  else
    error('Value must be either 0, 1, 2, NaN, or an array.');
  end
elseif isvector(x)
  if all(isfinite(x)) && all(x >= 0)
    % Input valid
  else
    error('Array must contain only finite, nonnegative numbers.');
  end
end
end

function mustBeNonnegativeOrNaN(x)
x = x(:);
if ~any(isnan(x)) && (~any(isfinite(x)) || any(x < 0))
  error('Value must be nonnegative or NaN.');
end
end