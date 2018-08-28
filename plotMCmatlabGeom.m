function plotMCmatlabGeom(G)
%   Created 2018 by Dominik Marti and Anders K. Hansen, DTU Fotonik
%
%   This function was inspired by lookmcxyz.m of the mcxyz MC program hosted at omlc.org
%
%   Input
%       name
%           the basename of the data files
%
%   Displays
%       Geometry cuboid
%       Media optical, thermal and fluorescence properties
%   If Monte Carlo output data exists, displays
%       Absorbed power
%       Fluence rate
%   And, if fluorescence Monte Carlo output data exists, displays
%       Media optical and thermal properties at the fluorescence wavelength
%       Distribution of fluorescence emitters
%       Absorbed fluorescence power
%       Fluorescence fluence rate
%	And, if light collectors were defined, displays
%		An illustration of the light collector angle and placement
%		Image generated
%
%   Requires
%       plotVolumetric.m
%       plotMediaProperties.m
%

%% Make geometry plot
h_f = plotVolumetric(1,G.x,G.y,G.z,G.M,'MCmatlab_GeometryIllustration',G.mediaProperties);
h_f.Name = 'Geometry illustration';
title('Geometry illustration');

%% Make media properties plot
h_f = plotMediaProperties(2,G.mediaProperties);
h_f.Name = 'Media properties';

if(~isnan(G.wavelength_f))
    %% Make fluorescence media properties plot
    h_f = plotMediaProperties(3,G.mediaProperties_f);
    h_f.Name = 'Fluorescence media properties';
end
drawnow;
end
