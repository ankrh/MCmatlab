addpath([fileparts(matlab.desktop.editor.getActiveFilename) '/helperfuncs']); % The helperfuncs folder is added to the path for the duration of this MATLAB session

%% Description
% Here, we demonstrate calculation and plotting of the far field of the
% "escaped" photons, both for excitation light and for fluorescence light.
% It is enabled by specifying the optional "farfieldRes" parameter, which
% serves also to specify the resolution you want of the far field
% distribution. The geometry is the same fluorescing cylinder as example 4,
% but now illuminated by a tilted Gaussian beam.
% 
% In the far field of the excitation light, you can see that they primarily
% escape in downward-pointing directions (in transmission), while the far
% field distribution of fluorescence light indicates that fluorescence
% mostly escapes on the ends of the cylinder, with a small amount of light
% coming out of the cylinder in the minus z direction (upwards).
% 
% Note that MCmatlab distinguishes between "escaped" photons and "killed"
% photons. An "escaping" photon is one that hits the top cuboid boundary
% (if boundaryType == 2) or any cuboid boundary (if boundaryType == 1)
% where the medium has refractive index 1 (or matchedInterfaces == true). A
% "killed" photon is one that strays too far from the main cuboid (6 times
% further than the cuboid dimensions).
% 
% If boundaryType == 1, there are no "killed" photons since no photons can
% travel outside the cuboid, and the fraction of light absorbed in the
% cuboid plus the fraction of light escaping equals 1.

%% Geometry definition
clear Ginput
Ginput.nx                = 101; % Number of bins in the x direction
Ginput.ny                = 101; % Number of bins in the y direction
Ginput.nz                = 101; % Number of bins in the z direction
Ginput.Lx                = .1; % [cm] x size of simulation cuboid
Ginput.Ly                = .1; % [cm] y size of simulation cuboid
Ginput.Lz                = .1; % [cm] z size of simulation cuboid

Ginput.mediaPropertiesFunc = @mediaPropertiesFunc; % Media properties defined as a function at the end of this file
Ginput.GeomFunc          = @GeometryDefinition_FluorescingCylinder; % Function to use for defining the distribution of media in the cuboid. Defined at the end of this m file.

% Execution, do not modify the next line:
Goutput = defineGeometry(Ginput);

plotMCmatlabGeom(Goutput);

%% Monte Carlo simulation
clear MCinput
MCinput.simulationTime           = .1; % [min] Time duration of the simulation
MCinput.farfieldRes              = 50; % (Default: 0) If nonzero, photons that "escape" will have their energies tracked in a 2D angle distribution (theta,phi) array with theta and phi resolutions equal to this number. An "escaping" photon is one that hits the top cuboid boundary (if boundaryType == 2) or any cuboid boundary (if boundaryType == 1) where the medium has refractive index 1.

MCinput.matchedInterfaces        = true; % Assumes all refractive indices are 1
MCinput.boundaryType             = 1; % 0: No escaping boundaries, 1: All cuboid boundaries are escaping, 2: Top cuboid boundary only is escaping
MCinput.wavelength               = 450; % [nm] Excitation wavelength, used for determination of optical properties for excitation light

MCinput.Beam.beamType            = 3; % 0: Pencil beam, 1: Isotropically emitting point source, 2: Infinite plane wave, 3: Gaussian focus, Gaussian far field beam, 4: Gaussian focus, top-hat far field beam, 5: Top-hat focus, Gaussian far field beam, 6: Top-hat focus, top-hat far field beam, 7: Laguerre-Gaussian LG01 beam
MCinput.Beam.xFocus              = 0; % [cm] x position of focus
MCinput.Beam.yFocus              = 0; % [cm] y position of focus
MCinput.Beam.zFocus              = 0.03; % [cm] z position of focus
MCinput.Beam.theta               = pi/6; % [rad] Polar angle of beam center axis
MCinput.Beam.phi                 = -pi/2; % [rad] Azimuthal angle of beam center axis
MCinput.Beam.waist               = 0.015; % [cm] Beam waist 1/e^2 radius
MCinput.Beam.divergence          = 15/180*pi; % [rad] Beam divergence 1/e^2 half-angle of beam (for a diffraction limited Gaussian beam, this is G.wavelength*1e-9/(pi*MCinput.Beam.waist*1e-2))

% Execution, do not modify the next two lines:
MCinput.G = Goutput;
MCoutput = runMonteCarlo(MCinput);

plotMCmatlab(MCinput,MCoutput);

%% Fluorescence Monte Carlo
clear FMCinput
FMCinput.simulationTime           = .1; % [min] Time duration of the simulation
FMCinput.farfieldRes              = 50; % (Default: 0) If nonzero, photons that "escape" will have their energies tracked in a 2D angle distribution (theta,phi) array with theta and phi resolutions equal to this number. An "escaping" photon is one that hits the top cuboid boundary (if boundaryType == 2) or any cuboid boundary (if boundaryType == 1) where the medium has refractive index 1.

FMCinput.matchedInterfaces        = true; % Assumes all refractive indices are 1
FMCinput.boundaryType             = 1; % 0: No escaping boundaries, 1: All cuboid boundaries are escaping, 2: Top cuboid boundary only is escaping
FMCinput.wavelength               = 550; % [nm] Fluorescence wavelength, used for determination of optical properties for fluorescence light

% Execution, do not modify the next three lines:
FMCinput.G = Goutput;
FMCinput.MCoutput = MCoutput;
FMCoutput = runMonteCarloFluorescence(FMCinput);

plotMCmatlabFluorescence(FMCinput,FMCoutput);

%% Post-processing

%% Geometry function(s)
% A geometry function takes as input X,Y,Z matrices as returned by the
% "ndgrid" MATLAB function as well as any parameters the user may have
% provided in the definition of Ginput. It returns the media matrix M,
% containing numerical values indicating the media type (as defined in
% getMediaProperties) at each voxel location.
function M = GeometryDefinition_FluorescingCylinder(X,Y,Z,parameters)
cylinderradius  = 0.0100;
M = 1*ones(size(X)); % fill background with fluorescence absorber
M(Y.^2 + (Z - 3*cylinderradius).^2 < cylinderradius^2) = 2; % fluorescer
end

%% Media Properties function
% The media properties function defines all the optical and thermal
% properties of the media involved by constructing and returning a
% "mediaProperties" struct with various fields. As its input, the function
% takes the wavelength as well as any other parameters you might specify
% above in the model file, for example parameters that you might loop over
% in a for loop.
function mediaProperties = mediaPropertiesFunc(wavelength,parameters)
j=1;
mediaProperties(j).name  = 'test fluorescence absorber';
if(wavelength<500)
  mediaProperties(j).mua = 1;
  mediaProperties(j).mus = 100;
  mediaProperties(j).g   = 0.9;
else
  mediaProperties(j).mua = 100;
  mediaProperties(j).mus = 100;
  mediaProperties(j).g   = 0.9;
end
mediaProperties(j).n   = 1.3;

j=2;
mediaProperties(j).name  = 'test fluorescer';
if(wavelength<500)
  mediaProperties(j).mua = 100;
  mediaProperties(j).mus = 100;
  mediaProperties(j).g   = 0.9;

  mediaProperties(j).Y   = 0.5;
else
  mediaProperties(j).mua = 1;
  mediaProperties(j).mus = 100;
  mediaProperties(j).g   = 0.9;
end
mediaProperties(j).n   = 1.3;

end
