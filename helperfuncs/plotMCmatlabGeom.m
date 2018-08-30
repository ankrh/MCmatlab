function plotMCmatlabGeom(G)
%   Created 2018 by Dominik Marti and Anders K. Hansen, DTU Fotonik
%
%   Displays
%       Geometry cuboid
%       Media optical, thermal and fluorescence properties
%
%   Requires
%       plotVolumetric.m
%       plotMediaProperties.m
%

%% Make fluorescence media properties plot
if(~isnan(G.wavelength_f))
    h_f = plotMediaProperties(3,G.mediaProperties_f);
    h_f.Name = 'Fluorescence media properties';
end

%% Make media properties plot
h_f = plotMediaProperties(2,G.mediaProperties);
h_f.Name = 'Media properties';

%% Make geometry plot
h_f = plotVolumetric(1,G.x,G.y,G.z,G.M,'MCmatlab_GeometryIllustration',G.mediaProperties);
h_f.Name = 'Geometry illustration';
title('Geometry illustration');

drawnow;
end
