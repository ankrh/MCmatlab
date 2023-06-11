%% Description
% This example shows how to import a shape from an STL file (sample.stl)
% and optionally apply scaling (unit conversion), rotation and mirroring on
% the mesh before using a function to voxelise the mesh.
%
% Importing is done inside the geometry function at the end of the model
% file (this file). The function findInsideVoxels returns a logical 3D
% array which is true in the indices of those voxels that are inside the
% (watertight) shape stored in the STL file. The function takes six inputs.
% The first three are always just the X, Y, and Z arrays exactly as they
% are passed into the geometry function. The fourth is the file name of the
% STL file. The fifth is a 3x3 matrix A that will be multiplied onto the
% STL mesh [x;y;z] coordinates. The sixth is a 1D array v of length 3 that
% defines the translation that will be applied to the mesh points after
% multiplication and before voxelisation. In other words, the mesh points
% are transformed according to    [x';y';z'] = A*[x;y;z] + v     and the
% voxelisation is then performed based on the mesh with positions
% [x';y';z'].
%
% The multiplication matrix A can contain a scaling factor, e.g., for
% converting from the units of the STL file to centimeters, which is the
% unit used in MCmatlab. It can also contain various rotation matrices, see
% https://en.wikipedia.org/wiki/Rotation_matrix#In_three_dimensions
% You may also mirror the coordinates by inverting the sign of the
% appropriate elements in the matrix.
%
% The example runs four separate sub-examples showing various
% transformations. Press enter to advance from one sub-example to the next.
% The MC simulation is not run in this example.
%
% Remember that since the z-axis in MCmatlab points down, the coordinate
% system is a left-handed coordinate system. Therefore a rotation that is
% clockwise in a right-handed coordinate system will be counter-clockwise
% in MCmatlab's coordinate system.
%
% The STL import feature uses the mesh_voxelization code by William Warriner:
% https://github.com/wwarriner/mesh_voxelization
% although the function has been lightly modified by Anders K. Hansen for
% use in MCmatlab.

%% MCmatlab abbreviations
% G: Geometry, MC: Monte Carlo, FMC: Fluorescence Monte Carlo, HS: Heat
% simulation, M: Media array, FR: Fluence rate, FD: Fractional damage.
%
% There are also some optional abbreviations you can use when referencing
% object/variable names: LS = lightSource, LC = lightCollector, FPID =
% focalPlaneIntensityDistribution, AID = angularIntensityDistribution, NI =
% normalizedIrradiance, NFR = normalizedFluenceRate.
%
% For example, "model.MC.LS.FPID.radialDistr" is the same as
% "model.MC.lightSource.focalPlaneIntensityDistribution.radialDistr"

%% Geometry definition
MCmatlab.closeMCmatlabFigures();
model = MCmatlab.model;

model.G.nx                = 100; % Number of bins in the x direction
model.G.ny                = 100; % Number of bins in the y direction
model.G.nz                = 100; % Number of bins in the z direction
model.G.Lx                = .3; % [cm] x size of simulation cuboid
model.G.Ly                = .3; % [cm] y size of simulation cuboid
model.G.Lz                = .3; % [cm] z size of simulation cuboid

model.G.mediaPropertiesFunc = @mediaPropertiesFunc; % Media properties defined as a function at the end of this file

model.G.geomFunc          = @geometryDefinition_1;
model = plot(model,'G');
% figure(31); % Make the STL import illustration the active figure so we can see it


%% Monte Carlo simulation
model.MC.simulationTimeRequested  = .1; % [min] Time duration of the simulation
model.MC.matchedInterfaces        = true; % Assumes all refractive indices are the same
model.MC.boundaryType             = 1; % 0: No escaping boundaries, 1: All cuboid boundaries are escaping, 2: Top cuboid boundary only is escaping, 3: Top and bottom boundaries are escaping, while the side boundaries are cyclic
model.MC.wavelength               = 532; % [nm] Excitation wavelength, used for determination of optical properties for excitation light

model.MC.lightSource.sourceType   = 4; % 0: Pencil beam, 1: Isotropically emitting line or point source, 2: Infinite plane wave, 3: Laguerre-Gaussian LG01 beam, 4: Radial-factorizable beam (e.g., a Gaussian beam), 5: X/Y factorizable beam (e.g., a rectangular LED emitter)
model.MC.lightSource.focalPlaneIntensityDistribution.radialDistr = 0; % Radial focal plane intensity distribution - 0: Top-hat, 1: Gaussian, Array: Custom. Doesn't need to be normalized.
model.MC.lightSource.focalPlaneIntensityDistribution.radialWidth = .03; % [cm] Radial focal plane 1/e^2 radius if top-hat or Gaussian or half-width of the full distribution if custom
model.MC.lightSource.angularIntensityDistribution.radialDistr = 0; % Radial angular intensity distribution - 0: Top-hat, 1: Gaussian, 2: Cosine (Lambertian), Array: Custom. Doesn't need to be normalized.
model.MC.lightSource.angularIntensityDistribution.radialWidth = 0; % [rad] Radial angular 1/e^2 half-angle if top-hat or Gaussian or half-angle of the full distribution if custom. For a diffraction limited Gaussian beam, this should be set to model.MC.wavelength*1e-9/(pi*model.MC.lightSource.focalPlaneIntensityDistribution.radialWidth*1e-2))
model.MC.lightSource.xFocus       = 0; % [cm] x position of focus
model.MC.lightSource.yFocus       = 0; % [cm] y position of focus
model.MC.lightSource.zFocus       = 0; % [cm] z position of focus
model.MC.lightSource.theta        = 0; % [rad] Polar angle of beam center axis
model.MC.lightSource.phi          = 0; % [rad] Azimuthal angle of beam center axis


% These lines will run the Monte Carlo simulation with the provided
% parameters and subsequently plot the results:
model = runMonteCarlo(model);
model = plot(model,'MC');




%% Geometry function(s) (see readme for details)
function M = geometryDefinition_1(X,Y,Z,parameters)
  theta = pi/2;
  A = 0.1*[cos(theta)   0   sin(theta)    ;
               0        1      0;
           -sin(theta)  0   cos(theta)]; % Rotation around the y axis

  % Rotate STL mesh points, convert the coordinates from mm to cm by multiplying by 0.1:
  A = 0.1*A;

  % Then translate the points by 1.5 in the +z direction:
  v = [-0.25 -0.4 0.55];

  % Find out which voxels are located inside this mesh:
  insideVoxels = findInsideVoxels(X,Y,Z,'brain.stl',A,v);

  % Set the background to air and the inside voxels to standard tissue:
  M = ones(size(X)); % Air
  M(insideVoxels) = 2;
end
%% Media Properties function (see readme for details)
function mediaProperties = mediaPropertiesFunc(parameters)
  mediaProperties = MCmatlab.mediumProperties;

  j=1;
  mediaProperties(j).name  = 'air';
  mediaProperties(j).mua   = 1e-8; % [cm^-1]
  mediaProperties(j).mus   = 1e-8; % [cm^-1]
  mediaProperties(j).g     = 1;

  j=2;
  mediaProperties(j).name  = 'brain tissue';
  mediaProperties(j).mua   = 0.1; % [cm^-1]
  mediaProperties(j).mus   = 200; % [cm^-1]
  mediaProperties(j).g     = 0.9;
end