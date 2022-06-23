%% Description
% In this introductory example, a block of "standard tissue" (mu_a = 1,
% mu_s = 100, g = 0.9) is illuminated by a pencil beam (infinitely thin
% beam). A small slice of air is present in the top of the simulation
% volume. nx and ny are set to odd values so that when the pencil beam is
% launched at x = y = 0 and travels straight down, it travels along a
% well-defined center column of voxels (the middle of the 51st column). Use
% the log10 plot checkbox in the visualizations to better see the fluence
% rate and absorption distribution in the MC result.

%% Geometry definition
model = MCmatlab.model;

model.G.nx                = 101; % Number of bins in the x direction
model.G.ny                = 101; % Number of bins in the y direction
model.G.nz                = 150; % Number of bins in the z direction
model.G.Lx                = .1; % [cm] x size of simulation cuboid
model.G.Ly                = .1; % [cm] y size of simulation cuboid
model.G.Lz                = .15; % [cm] z size of simulation cuboid

model.G.mediaPropertiesFunc = @mediaPropertiesFunc; % Media properties defined as a function at the end of this file
model.G.geomFunc          = @geometryDefinition; % Function to use for defining the distribution of media in the cuboid. Defined at the end of this m file.

plot(model,'G');

%% Monte Carlo simulation
model.MC.simulationTimeRequested  = .1; % [min] Time duration of the simulation
model.MC.matchedInterfaces        = true; % Assumes all refractive indices are the same
model.MC.boundaryType             = 1; % 0: No escaping boundaries, 1: All cuboid boundaries are escaping, 2: Top cuboid boundary only is escaping
model.MC.wavelength               = 532; % [nm] Excitation wavelength, used for determination of optical properties for excitation light

model.MC.lightSource.sourceType   = 0; % 0: Pencil beam, 1: Isotropically emitting line or point source, 2: Infinite plane wave, 3: Laguerre-Gaussian LG01 beam, 4: Radial-factorizable beam (e.g., a Gaussian beam), 5: X/Y factorizable beam (e.g., a rectangular LED emitter)
model.MC.lightSource.xFocus       = 0; % [cm] x position of focus
model.MC.lightSource.yFocus       = 0; % [cm] y position of focus
model.MC.lightSource.zFocus       = model.G.Lz/2; % [cm] z position of focus
model.MC.lightSource.theta        = 0; % [rad] Polar angle of beam center axis
model.MC.lightSource.phi          = 0; % [rad] Azimuthal angle of beam center axis


% Execution, do not modify the next line:
model = runMonteCarlo(model);

plot(model,'MC');

%% Geometry function(s)
% A geometry function takes as input X,Y,Z matrices as returned by the
% "ndgrid" MATLAB function as well as any parameters the user may have
% provided in the definition of Ginput. It returns the media matrix M,
% containing numerical values indicating the media type (as defined in
% mediaPropertiesFunc) at each voxel location.
function M = geometryDefinition(X,Y,Z,parameters)
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
% fractional heat damage FD can be specified as in examples 12-15.
function mediaProperties = mediaPropertiesFunc(wavelength,parameters)
    j=1;
    mediaProperties(j).name  = 'air';
    mediaProperties(j).mua   = 1e-8; % [cm^-1]
    mediaProperties(j).mus   = 1e-8; % [cm^-1]
    mediaProperties(j).g     = 1;
    mediaProperties(j).n     = 1;
    
    j=2;
    mediaProperties(j).name  = 'standard tissue';
    mediaProperties(j).mua   = 1; % [cm^-1]
    mediaProperties(j).mus   = 100; % [cm^-1]
    mediaProperties(j).g     = 0.9;
    mediaProperties(j).n     = 1.3;
end