function [M_muasplit, mediaProperties] = getMediaProperties(M_unsplit,mediaProperties_funcHandles,I,T)
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
n_muaBins = 5;
n_musBins = 6;
n_gBins = 7;

%% Split the mua dependent media into bins
M_muasplit = NaN(size(M_unsplit));
mediaProperties_muasplit = mediaProperties_funcHandles;
j = 1; % Position in mediaProperties_muasplit
for i=1:length(mediaProperties_funcHandles)
  if isa(mediaProperties_funcHandles(i).mua,'function_handle')
    mua_values = mediaProperties_funcHandles(i).mua(I(M_unsplit(:) == i),T(M_unsplit(:) == i));
    bin_values = linspace(min(mua_values),max(mua_values),n_muaBins);
    M_muasplit(M_unsplit(:) == i) = round((mua_values - min(mua_values))/(max(mua_values)-min(mua_values))*(n_muaBins-1)) + i;
    for k=0:n_muaBins-1
      mediaProperties_muasplit(j+k) = mediaProperties_funcHandles(i);
      mediaProperties_muasplit(j+k).mua = bin_values(k+1);
    end
    j = j + n_muaBins;
  else
    M_muasplit(M_unsplit(:) == i) = i;
    mediaProperties_muasplit(j) = mediaProperties_funcHandles(i);
    j = j + 1;
  end
end

%% Split the mus dependent media into bins
M_musmuasplit = NaN(size(M_muasplit));
mediaProperties_musmuasplit = mediaProperties_funcHandles;
j = 1; % Position in mediaProperties_musmuasplit
for i=1:length(mediaProperties_muasplit)
  if isa(mediaProperties_muasplit(i).mus,'function_handle')
    mus_values = mediaProperties_muasplit(i).mus(I(M_muasplit(:) == i),T(M_muasplit(:) == i));
    bin_values = linspace(min(mus_values),max(mus_values),n_musBins);
    M_musmuasplit(M_muasplit(:) == i) = round((mus_values - min(mus_values))/(max(mus_values)-min(mus_values))*(n_musBins-1)) + i;
    for k=0:n_musBins-1
      mediaProperties_musmuasplit(j+k) = mediaProperties_muasplit(i);
      mediaProperties_musmuasplit(j+k).mus = bin_values(k+1);
    end
    j = j + n_musBins;
  else
    M_musmuasplit(M_muasplit(:) == i) = i;
    mediaProperties_musmuasplit(j) = mediaProperties_muasplit(i);
    j = j + 1;
  end
end

%% Split the g dependent media into bins
M_gmusmuasplit = NaN(size(M_unsplit));
mediaProperties_gmusmuasplit = mediaProperties_musmuasplit;
j = 1; % Position in mediaProperties_gmusmuasplit
for i=1:length(mediaProperties_musmuasplit)
  if isa(mediaProperties_musmuasplit(i).g,'function_handle')
    g_values = mediaProperties_musmuasplit(i).g(I(M_musmuasplit(:) == i),T(M_musmuasplit(:) == i));
    bin_values = linspace(min(g_values),max(g_values),n_gBins);
    M_gmusmuasplit(M_musmuasplit(:) == i) = round((g_values - min(g_values))/(max(g_values)-min(g_values))*(n_gBins-1)) + i;
    for k=0:n_gBins-1
      mediaProperties_gmusmuasplit(j+k) = mediaProperties_musmuasplit(i);
      mediaProperties_gmusmuasplit(j+k).g = bin_values(k+1);
    end
    j = j + n_gBins;
  else
    M_gmusmuasplit(M_musmuasplit(:) == i) = i;
    mediaProperties_gmusmuasplit(j) = mediaProperties_musmuasplit(i);
    j = j + 1;
  end
end

uniqueParameterCombinations = unique([mua_mat(:) mus_mat(:) g_mat(:)],'rows');

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
  elseif(isfield(mediaProperties,'n') && ~isempty(mediaProperties(j).n) && mediaProperties(j).n < 1)
    error('Medium %s has n < 1',mediaProperties(j).name);
  elseif(isfield(mediaProperties,'VHC') && ~isempty(mediaProperties(j).VHC) && mediaProperties(j).VHC <= 0)
    error('Medium %s has VHC <= 0',mediaProperties(j).name);
  elseif(isfield(mediaProperties,'TC') && ~isempty(mediaProperties(j).TC) && mediaProperties(j).TC < 0)
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