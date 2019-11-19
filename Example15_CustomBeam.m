addpath([fileparts(matlab.desktop.editor.getActiveFilename) '/helperfuncs']); % The helperfuncs folder is added to the path for the duration of this MATLAB session
fprintf('\n');

%% Description
% 

%% Geometry definition
model = initializeMCmatlabModel();

model.G.nx                = 101; % Number of bins in the x direction
model.G.ny                = 101; % Number of bins in the y direction
model.G.nz                = 150; % Number of bins in the z direction
model.G.Lx                = .1; % [cm] x size of simulation cuboid
model.G.Ly                = .1; % [cm] y size of simulation cuboid
model.G.Lz                = .15; % [cm] z size of simulation cuboid

model.G.mediaPropertiesFunc = @mediaPropertiesFunc; % Media properties defined as a function at the end of this file
model.G.geomFunc          = @geometryDefinition_StandardTissue; % Function to use for defining the distribution of media in the cuboid. Defined at the end of this m file.

% Execution, do not modify the next line:
model = defineGeometry(model);

plotMCmatlabGeom(model);

%% Monte Carlo simulation
model = clearMCmatlabModel(model,'MC'); % Only necessary if you want to run this section repeatedly, re-using previous G data

model.MC.simulationTime           = .1; % [min] Time duration of the simulation

model.MC.matchedInterfaces        = true; % Assumes all refractive indices are 1
model.MC.boundaryType             = 1; % 0: No escaping boundaries, 1: All cuboid boundaries are escaping, 2: Top cuboid boundary only is escaping
model.MC.wavelength               = 532; % [nm] Excitation wavelength, used for determination of optical properties for excitation light

model.MC.beam.beamType            = 0; % 0: Pencil beam, 1: Isotropically emitting point source, 2: Infinite plane wave, 3: Laguerre-Gaussian LG01 beam, 4: Custom Radial/Azimuthal beam, 5: Custom X/Y beam

model.MC.beam.NF.radialDistr      = 1; % Radial near field distribution - 0: Top-hat, 1: Gaussian, Array: Custom. Doesn't need to be normalized.
model.MC.beam.NF.radialWidth      = .01; % [cm] Radial near field half-width of the full distribution if custom or 1/e^2 radius if top-hat or Gaussian
model.MC.beam.NF.azimuthalDistr   = 1; % Azimuthal near field distribution (which is assumed to span the angles (-pi, pi]) - 0: Uniform, Array: Custom. Doesn't need to be normalized.
model.MC.beam.NF.XDistr           = 1; % X near field distribution - 0: Top-hat, 1: Gaussian, Array: Custom. Doesn't need to be normalized.
model.MC.beam.NF.XWidth           = .01; % [cm] X near field half-width of the full distribution if custom or 1/e^2 radius if top-hat or Gaussian
model.MC.beam.NF.YDistr           = 1; % Y near field distribution - 0: Top-hat, 1: Gaussian, Array: Custom. Doesn't need to be normalized.
model.MC.beam.NF.YWidth           = .01; % [cm] Y near field half-width of the full distribution if custom or 1/e^2 radius if top-hat or Gaussian

model.MC.beam.FF.radialDistr      = 1; % Radial far field distribution - 0: Top-hat, 1: Gaussian, 2: Cosine (Lambertian), Array: Custom. Doesn't need to be normalized.
model.MC.beam.FF.radialWidth      = pi/4; % [rad] Radial far field half-angle of the full distribution if custom or 1/e^2 half-angle if top-hat or Gaussian. For a diffraction limited Gaussian beam, this should be set to model.MC.wavelength*1e-9/(pi*model.MC.beam.NF.radialWidth*1e-2))
model.MC.beam.FF.azimuthalDistr   = 1; % Azimuthal far field distribution (which is assumed to span the angles (-pi, pi]) - 0: Uniform, Array: Custom. Doesn't need to be normalized.
model.MC.beam.FF.XDistr           = 1; % X far field distribution - 0: Top-hat, 1: Gaussian, Array: Custom. Doesn't need to be normalized.
model.MC.beam.FF.XWidth           = pi/4; % [rad] X far field half-angle of the full distribution if custom or 1/e^2 half-angle if top-hat or Gaussian
model.MC.beam.FF.YDistr           = 1; % Y far field distribution - 0: Top-hat, 1: Gaussian, Array: Custom. Doesn't need to be normalized.
model.MC.beam.FF.YWidth           = pi/4; % [rad] Y far field half-angle of the full distribution if custom or 1/e^2 half-angle if top-hat or Gaussian

model.MC.beam.xFocus              = 0; % [cm] x position of focus
model.MC.beam.yFocus              = 0; % [cm] y position of focus
model.MC.beam.zFocus              = model.G.Lz/2; % [cm] z position of focus

model.MC.beam.theta               = 0; % [rad] Polar angle of beam center axis
model.MC.beam.phi                 = 0; % [rad] Azimuthal angle of beam center axis
model.MC.beam.psi                 = 0; % [rad] Axial rotation angle of beam, relevant only for XY distributed beams or beams with nonuniform azimuthal distribution

% Execution, do not modify the next line:
model = runMonteCarlo(model);

plotMCmatlab(model);

%% Post-processing

%% Geometry function(s)
% A geometry function takes as input X,Y,Z matrices as returned by the
% "ndgrid" MATLAB function as well as any parameters the user may have
% provided in the definition of Ginput. It returns the media matrix M,
% containing numerical values indicating the media type (as defined in
% mediaPropertiesFunc) at each voxel location.
function M = geometryDefinition_StandardTissue(X,Y,Z,parameters)
tissuedepth = 0.03;
M = ones(size(X)); % Air
M(Z > tissuedepth) = 2; % "Standard" tissue
end

%% Media Properties function
% The media properties function defines all the optical and thermal
% properties of the media involved by constructing and returning a
% "mediaProperties" struct with various fields. As its input, the function
% takes the wavelength as well as any other parameters you might specify
% above in the model file, for example parameters that you might loop over
% in a for loop. Dependence on excitation fluence rate FR, temperature T or
% fractional heat damage FD can be specified as in examples 11-14.
function mediaProperties = mediaPropertiesFunc(wavelength,parameters)
j=1;
mediaProperties(j).name  = 'air';
mediaProperties(j).mua   = 1e-8;
mediaProperties(j).mus   = 1e-8;
mediaProperties(j).g     = 1;
mediaProperties(j).n     = 1;
mediaProperties(j).VHC   = 1.2e-3;
mediaProperties(j).TC    = 0; % Real value is 2.6e-4, but we set it to zero to neglect the heat transport to air

j=2;
mediaProperties(j).name  = 'standard tissue';
mediaProperties(j).mua   = 1;
mediaProperties(j).mus   = 100;
mediaProperties(j).g     = 0.9;
mediaProperties(j).n     = 1.3;
mediaProperties(j).VHC   = 3391*1.109e-3;
mediaProperties(j).TC    = 0.37e-2;
end