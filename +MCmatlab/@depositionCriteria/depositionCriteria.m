classdef depositionCriteria
% List of criteria for when a photon's weight loss should be deposited in
% any output arrays. All the criteria have to be met simultaneously - there
% is an implied AND between them all.
  properties
    minScatterings (1,1) double {mustBeInteger,mustBeNonnegative} = 0
    maxScatterings (1,1) double {mustBeNonnegativeIntegerOrPositiveInfinite} = Inf
    minRefractions (1,1) double {mustBeInteger,mustBeNonnegative} = 0
    maxRefractions (1,1) double {mustBeNonnegativeIntegerOrPositiveInfinite} = Inf
    minReflections (1,1) double {mustBeInteger,mustBeNonnegative} = 0
    maxReflections (1,1) double {mustBeNonnegativeIntegerOrPositiveInfinite} = Inf
    minInterfaceTransitions (1,1) double {mustBeInteger,mustBeNonnegative} = 0
    maxInterfaceTransitions (1,1) double {mustBeNonnegativeIntegerOrPositiveInfinite} = Inf
  end
end

function mustBeNonnegativeIntegerOrPositiveInfinite(x)
  if x < 0 || (isfinite(x) && rem(x,1) ~= 0)
    error('Error: Value must be a non-negative integer or infinity.');
  end
end