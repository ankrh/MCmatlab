%% Description
% In this example, we show two things: The use of cyclic boundary
% conditions and how to calculate the reflectance of a sample, including
% how to provide an estimate of the statistical error of the value.
%
% The geometry is similar to that of example 1, although we will run it
% both with and without matched interfaces here. We place a thin layer of
% air, 1 voxel thick, at the top of the cuboid and fill the rest with
% "standard tissue". The thin layer of air is necessary to enforce specular
% reflectance of the incident beam and to include the reflection and
% refraction effects at the surface.
%
% We want to calculate the total and diffuse reflectance of the tissue, so
% we don't want photons to escape at the side walls and disappear from the
% simulation. One possibility for how to avoid this is to use a cuboid with
% very large Lx and Ly. Another is to use cyclic boundary conditions, which
% we emply here by setting model.MC.boundaryType = 3. With cyclic boundary
% conditions, all photons that hit a side wall will immediately enter the
% cuboid again on the opposite wall. Because our geometry is supposed to
% represent a slab of tissue with infinite horizontal extent, this is a
% valid way for us to avoid losing photons to side wall effects.
%
% We will calculate the reflectance by integrating the
% model.MC.normalizedIrradiance_zneg 2D array, which contains all the power
% that has hit the top cuboid boundary from the inside. In principle, the x
% and y resolution can be set as low as the minimum value, 2x2 and still
% get the same result (try it). For visualization purposes, however, we
% keep nx and ny at reasonable values of 101 in this example. Also Lx and
% Ly could be set arbitrarily low due to the use of the cyclic boundaryType.
%
% The reflectance that we calculate is the total refletance, including the
% specular reflection that we get when simulating without matched
% interfaces. We can obtain the diffuse reflectance by subtracting the
% Fresnel reflectivity, ((1 - 1.4)/(1 + 1.4))^2, from the total
% reflectance. Another way of avoiding the specular reflectance would
% have been to set depositionCriteria.minInterfaceTransitions = 1, but we
% do not use that method here.
%
% To get some statistics, we launch 1e6 photons 5 times and collect the
% diffuse reflectance value for each run. Then we finally write out the
% mean and the standard error of the mean for both the calculation with
% matched interfaces and without.

%% Geometry definition
model = MCmatlab.model;

model.G.nx                = 101; % Number of bins in the x direction
model.G.ny                = 101; % Number of bins in the y direction
model.G.nz                = 100; % Number of bins in the z direction
model.G.Lx                = .5; % [cm] x size of simulation cuboid
model.G.Ly                = .5; % [cm] y size of simulation cuboid
model.G.Lz                = 2; % [cm] z size of simulation cuboid

model.G.mediaPropertiesFunc = @mediaPropertiesFunc; % Media properties defined as a function at the end of this file
model.G.geomFunc          = @geometryDefinition; % Function to use for defining the distribution of media in the cuboid. Defined at the end of this m file.

model = plot(model,'G');

%% Monte Carlo simulation
model.MC.nPhotonsRequested        = 1e5; % Number of photons requested
model.MC.boundaryType             = 3; % 0: No escaping boundaries, 1: All cuboid boundaries are escaping, 2: Top cuboid boundary only is escaping
model.MC.wavelength               = 532; % [nm] Excitation wavelength, used for determination of optical properties for excitation light

model.MC.lightSource.sourceType          = 0; % 0: Pencil beam, 1: Isotropically emitting line or point source, 2: Infinite plane wave, 3: Laguerre-Gaussian LG01 beam, 4: Radial-factorizable beam (e.g., a Gaussian beam), 5: X/Y factorizable beam (e.g., a rectangular LED emitter)
model.MC.lightSource.xFocus              = 0; % [cm] x position of focus
model.MC.lightSource.yFocus              = 0; % [cm] y position of focus
model.MC.lightSource.zFocus              = 0; % [cm] z position of focus

% Set up arrays to collect statistics in:
n = 5; % Times to run
Rd_matched = NaN(1,n);
Rd_mismatched = NaN(1,n);

% First run with matched interfaces:
model.MC.matchedInterfaces = true; % Assumes all refractive indices are the same
tic
for iRun = 1:n
  model = runMonteCarlo(model);

  Rd_matched(iRun) = sum(model.MC.normalizedIrradiance_zneg(:))*model.G.dx*model.G.dy; % Diffuse reflectance
end
t_Matched = toc;

% Then run with mismatched interfaces:
model.MC.matchedInterfaces = false;
tic
for iRun = 1:n
  model = runMonteCarlo(model);

  Rt_mismatched = sum(model.MC.normalizedIrradiance_zneg(:))*model.G.dx*model.G.dy; % Total reflectance
  Rd_mismatched(iRun) = Rt_mismatched - ((1 - 1.4)/(1 + 1.4))^2; % Diffuse reflectance
end
t_Mismatched = toc;

model = plot(model,'MC');

%% Print outputs:
fprintf('\nRd matched    = %.6f +- %.6f (n = %d, total time elapsed = %d s)\n'   ,mean(Rd_matched   ),std(Rd_matched   )/sqrt(n),n,round(t_Matched   ));
fprintf('Rd mismatched = %.6f +- %.6f (n = %d, total time elapsed = %d s)\n\n',mean(Rd_mismatched),std(Rd_mismatched)/sqrt(n),n,round(t_Mismatched));

%% Geometry function(s)
% A geometry function takes as input X,Y,Z matrices as returned by the
% "ndgrid" MATLAB function as well as any parameters the user may have
% provided in the definition of Ginput. It returns the media matrix M,
% containing numerical values indicating the media type (as defined in
% mediaPropertiesFunc) at each voxel location.
function M = geometryDefinition(X,Y,Z,parameters)
  M = 2*ones(size(X)); % "Standard" tissue
  M(:,:,1) = 1; % air
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
  mediaProperties(j).mus   = 100; % [cm^-1]
  mediaProperties(j).g     = 1;
  mediaProperties(j).n     = 1;

  j=2;
  mediaProperties(j).name  = 'standard tissue';
  mediaProperties(j).mua   = 1; % [cm^-1]
  mediaProperties(j).mus   = 100; % [cm^-1]
  mediaProperties(j).g     = 0.90;
  mediaProperties(j).n     = 1.4;
end