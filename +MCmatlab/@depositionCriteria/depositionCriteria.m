classdef depositionCriteria
% List of criteria for when a photon's weight loss should be deposited in
% any output arrays. All the criteria have to be met simultaneously - there
% is an implied AND between them all.

% The minMediumIdxToConsider and maxMediumIdxToConsider are different in
% that they specify which media the above criteria are for. For example, if
% you've specified minScatterings = 2, minMediumIdxToConsider = 4 and
% maxMediumIdxToConsider = 5, then only photons that have experienced at
% least 2 scattering events in media 4 and/or 5 will deposit their
% weight in the output arrays. It doesn't matter how many scattering events
% have happened in other media.

% Interface transitions and refractions will be considered if the medium
% transitioned *into* has index between minMediumIdxToConsider and
% maxMediumIdxToConsider.

% The user can rearrange the media in the media definition as required so
% that the set of media that he/she wants to consider in deposition
% criteria are in a contiguous interval.
  properties
    minScatterings (1,1) double {mustBeInteger,mustBeNonnegative} = 0
    maxScatterings (1,1) double {mustBeNonnegativeIntegerOrPositiveInfinite} = Inf
    minRefractions (1,1) double {mustBeInteger,mustBeNonnegative} = 0
    maxRefractions (1,1) double {mustBeNonnegativeIntegerOrPositiveInfinite} = Inf
    minReflections (1,1) double {mustBeInteger,mustBeNonnegative} = 0
    maxReflections (1,1) double {mustBeNonnegativeIntegerOrPositiveInfinite} = Inf
    minInterfaceTransitions (1,1) double {mustBeInteger,mustBeNonnegative} = 0
    maxInterfaceTransitions (1,1) double {mustBeNonnegativeIntegerOrPositiveInfinite} = Inf

    minMediumIdxToConsider (1,1) double {mustBePositiveIntegerOrPositiveInfinite} = 1
    maxMediumIdxToConsider (1,1) double {mustBePositiveIntegerOrPositiveInfinite} = Inf
  end

  properties (Hidden) % Calculated based on minMediumIdxToConsider and maxMediumIdxToConsider
    minSubmediaIdx (1,1) double
    maxSubmediaIdx (1,1) double 
  end
end

function mustBePositiveIntegerOrPositiveInfinite(x)
  if x <= 0 || (isfinite(x) && rem(x,1) ~= 0)
    error('Error: Value must be a positive integer or infinity.');
  end
end
function mustBeNonnegativeIntegerOrPositiveInfinite(x)
  if x < 0 || (isfinite(x) && rem(x,1) ~= 0)
    error('Error: Value must be a non-negative integer or infinity.');
  end
end