%% Description
% Here we demonstrate how to define beams using beamType = 1, 4 or 5. The
% geometry is only air, which has negligible scattering and absorption. The
% farFieldRes input will be further illustrated in example 7.
% 
% This example will execute five different simulations. Press any key with
% the command window in focus to start the next simulation.
% 
% First we demonstrate an isotropically emitting line source with beamType
% = 1 with length specified using beam.emitterLength and focus in the
% middle of the cuboid. If emitterLength is set to 0, the source would be a
% point. The theta and phi parameters determine what direction the line is
% pointing.
%
% Next, we set the focus at x=y=z=0 and theta = 0 so that the beam goes
% straight down. Then we define various simple and complicated beams using
% beamType 4 and 5. See other examples for more beams with focus placed
% inside the cuboid and for tilted input beams.
% 
% In principle, Gaussian beams can be simulated using either beamType 4 or
% 5 (because exp(-R^2) = exp(-X^2)*exp(-Y^2)).

%% Geometry definition
model = MCmatlab.model;

model.G.nx                = 201; % Number of bins in the x direction
model.G.ny                = 201; % Number of bins in the y direction
model.G.nz                = 100; % Number of bins in the z direction
model.G.Lx                = .2; % [cm] x size of simulation cuboid
model.G.Ly                = .2; % [cm] y size of simulation cuboid
model.G.Lz                = .1; % [cm] z size of simulation cuboid

model.G.mediaPropertiesFunc = @mediaPropertiesFunc; % Media properties defined as a function at the end of this file
model.G.geomFunc          = @geometryDefinition_Air; % Function to use for defining the distribution of media in the cuboid. Defined at the end of this m file.

% plot(model, 'G');

%% Monte Carlo simulations
%% Isotropic line emitter
model.MC.simulationTimeRequested  = .1; % [min] Time duration of the simulation
model.MC.farFieldRes              = 200; % (Default: 0) If nonzero, photons that "escape" will have their energies tracked in a 2D angle distribution (theta,phi) array with theta and phi resolutions equal to this number. An "escaping" photon is one that hits the top cuboid boundary (if boundaryType == 2) or any cuboid boundary (if boundaryType == 1) where the medium has refractive index 1.

model.MC.matchedInterfaces        = true; % Assumes all refractive indices are the same
model.MC.boundaryType             = 1; % 0: No escaping boundaries, 1: All cuboid boundaries are escaping, 2: Top cuboid boundary only is escaping
model.MC.wavelength               = 532; % [nm] Excitation wavelength, used for determination of optical properties for excitation light

model.MC.beam.beamType            = 1; % 0: Pencil beam, 1: Isotropically emitting line or point source, 2: Infinite plane wave, 3: Laguerre-Gaussian LG01 beam, 4: Radial-factorizable beam (e.g., a Gaussian beam), 5: X/Y factorizable beam (e.g., a rectangular LED emitter)
model.MC.beam.emitterLength       = 0.05; % [cm] (Optional) Length of the isotropically emitting line source (if 0, the source is a point source)

model.MC.beam.xFocus              = 0; % [cm] x position of focus
model.MC.beam.yFocus              = 0; % [cm] y position of focus
model.MC.beam.zFocus              = 0.05; % [cm] z position of focus

model.MC.beam.theta               = pi/4; % [rad] Polar angle of beam center axis
model.MC.beam.phi                 = pi/2; % [rad] Azimuthal angle of beam center axis
model.MC.beam.psi                 = 0; % [rad] (Default: 0) Axial rotation angle of beam, relevant only for XY distributed beams

% Execution, do not modify the next line:
model = runMonteCarlo(model);
plot(model,'MC');
fprintf('Press enter to continue...\n');pause;

%% Top-hat near field, Gaussian far field
model.MC.beam.xFocus              = 0; % [cm] x position of focus
model.MC.beam.yFocus              = 0; % [cm] y position of focus
model.MC.beam.zFocus              = 0; % [cm] z position of focus

model.MC.beam.beamType            = 4; % 0: Pencil beam, 1: Isotropically emitting line or point source, 2: Infinite plane wave, 3: Laguerre-Gaussian LG01 beam, 4: Radial-factorizable beam (e.g., a Gaussian beam), 5: X/Y factorizable beam (e.g., a rectangular LED emitter)
model.MC.beam.NF.radialDistr      = 0; % Radial near field distribution - 0: Top-hat, 1: Gaussian, Array: Custom. Doesn't need to be normalized.
model.MC.beam.NF.radialWidth      = .01; % [cm] Radial near field 1/e^2 radius if top-hat or Gaussian or half-width of the full distribution if custom
model.MC.beam.FF.radialDistr      = 1; % Radial far field distribution - 0: Top-hat, 1: Gaussian, 2: Cosine (Lambertian), Array: Custom. Doesn't need to be normalized.
model.MC.beam.FF.radialWidth      = pi/6; % [rad] Radial far field 1/e^2 half-angle if top-hat or Gaussian or half-angle of the full distribution if custom. For a diffraction limited Gaussian beam, this should be set to model.MC.wavelength*1e-9/(pi*model.MC.beam.NF.radialWidth*1e-2))
model.MC.beam.theta               = 0; % [rad] Polar angle of beam center axis

% Execution, do not modify the next line:
model = runMonteCarlo(model);
plot(model,'MC');
fprintf('Press enter to continue...\n');pause;

%% LED-like emitter: Rectangular near field, Lambertian far field
model.MC.beam.beamType            = 5; % 0: Pencil beam, 1: Isotropically emitting line or point source, 2: Infinite plane wave, 3: Laguerre-Gaussian LG01 beam, 4: Radial-factorizable beam (e.g., a Gaussian beam), 5: X/Y factorizable beam (e.g., a rectangular LED emitter)
model.MC.beam.NF.XDistr           = 0; % X near field distribution - 0: Top-hat, 1: Gaussian, Array: Custom. Doesn't need to be normalized.
model.MC.beam.NF.XWidth           = .02; % [cm] X near field 1/e^2 radius if top-hat or Gaussian or half-width of the full distribution if custom
model.MC.beam.NF.YDistr           = 0; % Y near field distribution - 0: Top-hat, 1: Gaussian, Array: Custom. Doesn't need to be normalized.
model.MC.beam.NF.YWidth           = .01; % [cm] Y near field 1/e^2 radius if top-hat or Gaussian or half-width of the full distribution if custom
model.MC.beam.FF.XDistr           = 2; % X far field distribution - 0: Top-hat, 1: Gaussian, 2: Cosine (Lambertian), Array: Custom. Doesn't need to be normalized.
model.MC.beam.FF.XWidth           = pi/8; % [rad] X far field 1/e^2 half-angle if top-hat or Gaussian or half-angle of the full distribution if custom
model.MC.beam.FF.YDistr           = 2; % Y far field distribution - 0: Top-hat, 1: Gaussian, 2: Cosine (Lambertian), Array: Custom. Doesn't need to be normalized.
model.MC.beam.FF.YWidth           = pi/8; % [rad] Y far field 1/e^2 half-angle if top-hat or Gaussian or half-angle of the full distribution if custom
model.MC.beam.psi                 = pi/4; % [rad] (Default: 0) Axial rotation angle of beam, relevant only for XY distributed beams

% Execution, do not modify the next line:
model = runMonteCarlo(model);
plot(model,'MC');
fprintf('Press enter to continue...\n');pause;

%% Custom radial near field and far field
model.MC.beam.beamType            = 4; % 0: Pencil beam, 1: Isotropically emitting line or point source, 2: Infinite plane wave, 3: Laguerre-Gaussian LG01 beam, 4: Radial-factorizable beam (e.g., a Gaussian beam), 5: X/Y factorizable beam (e.g., a rectangular LED emitter)
model.MC.beam.NF.radialDistr      = exp(-(linspace(0,20,1000)-4).^2); % Radial near field distribution - 0: Top-hat, 1: Gaussian, Array: Custom. Doesn't need to be normalized.
model.MC.beam.NF.radialWidth      = .025; % [cm] Radial near field 1/e^2 radius if top-hat or Gaussian or half-width of the full distribution if custom
model.MC.beam.FF.radialDistr      = 1+cos(linspace(0,7*pi,1000)); % Radial far field distribution - 0: Top-hat, 1: Gaussian, 2: Cosine (Lambertian), Array: Custom. Doesn't need to be normalized.
model.MC.beam.FF.radialWidth      = pi/4; % [rad] Radial far field 1/e^2 half-angle if top-hat or Gaussian or half-angle of the full distribution if custom. For a diffraction limited Gaussian beam, this should be set to model.MC.wavelength*1e-9/(pi*model.MC.beam.NF.radialWidth*1e-2))

% Execution, do not modify the next line:
model = runMonteCarlo(model);
plot(model,'MC');
fprintf('Press enter to continue...\n');pause;

%% Custom XY near field, top-hat X far field and Gaussian Y far field
model.MC.beam.beamType            = 5; % 0: Pencil beam, 1: Isotropically emitting line or point source, 2: Infinite plane wave, 3: Laguerre-Gaussian LG01 beam, 4: Radial-factorizable beam (e.g., a Gaussian beam), 5: X/Y factorizable beam (e.g., a rectangular LED emitter)
model.MC.beam.NF.XDistr           = sin(linspace(0,2*pi,1000)).^2; % X near field distribution - 0: Top-hat, 1: Gaussian, Array: Custom. Doesn't need to be normalized.
model.MC.beam.NF.XWidth           = .02; % [cm] X near field 1/e^2 radius if top-hat or Gaussian or half-width of the full distribution if custom
model.MC.beam.NF.YDistr           = sin(linspace(0,3*pi,1000)).^2; % Y near field distribution - 0: Top-hat, 1: Gaussian, Array: Custom. Doesn't need to be normalized.
model.MC.beam.NF.YWidth           = .01; % [cm] Y near field 1/e^2 radius if top-hat or Gaussian or half-width of the full distribution if custom
model.MC.beam.FF.XDistr           = 0; % X far field distribution - 0: Top-hat, 1: Gaussian, 2: Cosine (Lambertian), Array: Custom. Doesn't need to be normalized.
model.MC.beam.FF.XWidth           = pi/8; % [rad] X far field 1/e^2 half-angle if top-hat or Gaussian or half-angle of the full distribution if custom
model.MC.beam.FF.YDistr           = 1; % Y far field distribution - 0: Top-hat, 1: Gaussian, 2: Cosine (Lambertian), Array: Custom. Doesn't need to be normalized.
model.MC.beam.FF.YWidth           = pi/8; % [rad] Y far field 1/e^2 half-angle if top-hat or Gaussian or half-angle of the full distribution if custom
model.MC.beam.psi                 = -pi/4; % [rad] (Default: 0) Axial rotation angle of beam, relevant only for XY distributed beams

% Execution, do not modify the next line:
model = runMonteCarlo(model);
plot(model,'MC');

%% Geometry function(s)
% A geometry function takes as input X,Y,Z matrices as returned by the
% "ndgrid" MATLAB function as well as any parameters the user may have
% provided in the definition of Ginput. It returns the media matrix M,
% containing numerical values indicating the media type (as defined in
% mediaPropertiesFunc) at each voxel location.
function M = geometryDefinition_Air(X,Y,Z,parameters)
    M = ones(size(X)); % Air
end

%% Media Properties function
% The media properties function defines all the optical and thermal
% properties of the media involved by constructing and returning a
% "mediaProperties" struct with various fields. As its input, the function
% takes the wavelength as well as any other parameters you might specify
% above in the model file, for example parameters that you might loop over
% in a for loop. Dependence on excitation fluence rate FR, temperature T or
% fractional heat damage FD can be specified as in examples 12-15.
function mediaProperties = mediaPropertiesFunc(wavelength,parameters)
    j=1;
    mediaProperties(j).name  = 'air';
    mediaProperties(j).mua   = 1e-8; % [cm^-1]
    mediaProperties(j).mus   = 1e-8; % [cm^-1]
    mediaProperties(j).g     = 1;
    mediaProperties(j).n     = 1;
    mediaProperties(j).VHC   = 1.2e-3;
    mediaProperties(j).TC    = 0; % Real value is 2.6e-4, but we set it to zero to neglect the heat transport to air
end