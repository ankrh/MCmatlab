function [M, mediaProperties] = getMediaProperties(M,wavelength,parameters)
%   Returns the reduced medium matrix, using only numbers from 1 up to the number of used media, and
%   the known media properties (optical, thermal and/or fluorescence) at the specified wavelength.
%
%   See also mediaPropertiesLibrary, defineGeometry

%%%%%
%   Copyright 2017, 2018 by Dominik Marti and Anders K. Hansen, DTU Fotonik
%   This function was inspired by makeTissueList.m of the mcxyz program hosted at omlc.org
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

mediaProperties = mediaPropertiesLibrary(wavelength,parameters);

%% Trim mediaProperties down to use only the media included in the input matrix M, and reduce M accordingly
mediumMap = zeros(1,length(mediaProperties),'uint8');
mediumMap(unique(M)) = 1:length(unique(M));
mediaProperties = mediaProperties(unique(M)); % Reduced medium list, containing only the used media
M = mediumMap(M); % Reduced medium matrix, using only numbers from 1 up to the number of used media

%% Fill in fluorescence and Arrhenius parameter assumptions
% For all media for which the fluorescence power yield Y, Arrhenius
% activation energy E or Arrhenius pre-exponential factor A was not 
% specified, assume they are zero.
for j=1:length(mediaProperties)
  if(~isfield(mediaProperties,'Y') || isempty(mediaProperties(j).Y))
    mediaProperties(j).Y = 0;
  end
  if(~isfield(mediaProperties,'E') || isempty(mediaProperties(j).E))
    mediaProperties(j).E = 0;
  end
  if(~isfield(mediaProperties,'A') || isempty(mediaProperties(j).A))
    mediaProperties(j).A = 0;
  end
end

%% Throw an error if a variable doesn't conform to its required interval
for j=1:length(mediaProperties)
  if(~isfield(mediaProperties,'mua') || isempty(mediaProperties(j).mua))
    error('Medium %s has no mua.',mediaProperties(j).name);
  elseif(~isfield(mediaProperties,'mus') || isempty(mediaProperties(j).mus))
    error('Medium %s has no mus.',mediaProperties(j).name);
  elseif(~isfield(mediaProperties,'g') || isempty(mediaProperties(j).g))
    error('Medium %s has no g.',mediaProperties(j).name);
  end

  if(mediaProperties(j).mua <= 0)
    error('Medium %s has mua <= 0',mediaProperties(j).name);
  elseif(mediaProperties(j).mus <= 0)
    error('Medium %s has mus <= 0',mediaProperties(j).name);
  elseif(abs(mediaProperties(j).g) > 1)
    error('Medium %s has abs(g) > 1',mediaProperties(j).name);
  elseif(mediaProperties(j).n < 1)
    error('Medium %s has n < 1',mediaProperties(j).name);
  elseif(mediaProperties(j).VHC <= 0)
    error('Medium %s has VHC <= 0',mediaProperties(j).name);
  elseif(mediaProperties(j).TC < 0)
    error('Medium %s has TC < 0',mediaProperties(j).name);
  elseif(mediaProperties(j).Y < 0)
    error('Medium %s has Y < 0',mediaProperties(j).name);
  elseif(mediaProperties(j).E < 0)
    error('Medium %s has E < 0',mediaProperties(j).name);
  elseif(mediaProperties(j).A < 0)
    error('Medium %s has A < 0',mediaProperties(j).name);
  end
end

end