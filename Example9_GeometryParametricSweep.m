%% Description
% This example shows how to execute MC simulations in a for loop, in this
% case simulating a pencil beam incident on a slab of "standard tissue"
% with variable (parametrically sweeped) thickness. The slab thickness is
% passed in as a part of the GeomFuncParams field and used in the geometry
% function to get the correct thickness. Light is collected in transmission
% at a 45° angle in a fiber. A small layer of air is present underneath the
% slab so that the refraction out of the slab is correctly modeled. At the
% end of the script, collected power as a function of thickness is plotted.
% The fiber-coupled power is seen to be zero for zero slab thickness, since
% there is nothing to scatter the light over into the fiber, and the power
% starts to drop off when the slab thickness passes 0.05 cm because then
% much of the light is either absorbed or scattered backwards rather than
% into the fiber.
%
% Lz and nz have been carefully chosen so that the slab interfaces always
% coincide with voxel boundaries, so we get exactly correct slab
% thicknesses in all the iterations of the for loop. Otherwise, the
% simulated slab thickness would deviate from the specified slab thickness
% because of the voxel rounding during the geometry definition.
%
% The silentMode flags are used here, which suppress the outputs to the
% command line, which is especially useful to avoid excessive text if
% simulating in a for- or while-loop like this.
% 
% Also, the optional calcNFR flag in the MCinput is here set to false, which
% means the MC simulation does not calculate the 3D fluence rate matrix.
% This is useful because we're here only interested in the "image" data,
% and setting calcF to false will speed up the simulation a bit (10-30%).


%% Geometry definition
model = MCmatlab.model;

model.G.silentMode          = true; % Disables command window text and progress indication

model.G.nx                  = 11; % Number of bins in the x direction
model.G.ny                  = 21; % Number of bins in the y direction
model.G.nz                  = 21; % Number of bins in the z direction
model.G.Lx                  = .1; % [cm] x size of simulation cuboid
model.G.Ly                  = .2; % [cm] y size of simulation cuboid
model.G.Lz                  = .105; % [cm] z size of simulation cuboid

model.G.mediaPropertiesFunc = @mediaPropertiesFunc; % Media properties defined as a function at the end of this file
model.G.geomFunc            = @geometryDefinition_VariableThicknessStandardTissue; % Function to use for defining the distribution of media in the cuboid. Defined at the end of this m file.

%% Monte Carlo simulation
model.MC.silentMode               = true; % Disables command window text and progress indication
model.MC.useAllCPUs               = true; % If false, MCmatlab will leave one processor unused. Useful for doing other work on the PC while simulations are running.
model.MC.simulationTimeRequested  = 2/60; % [min] Time duration of the simulation
model.MC.calcNFR                  = false; % (Default: true) If true, the 3D fluence rate output array NFR will be calculated. Set to false if you have a light collector and you're only interested in the image output.

model.MC.matchedInterfaces        = false; % Assumes all refractive indices are the same
model.MC.boundaryType             = 1; % 0: No escaping boundaries, 1: All cuboid boundaries are escaping, 2: Top cuboid boundary only is escaping
model.MC.wavelength               = 532; % [nm] Excitation wavelength, used for determination of optical properties for excitation light

model.MC.beam.beamType            = 0; % 0: Pencil beam, 1: Isotropically emitting line or point source, 2: Infinite plane wave, 3: Laguerre-Gaussian LG01 beam, 4: Radial-factorizable beam (e.g., a Gaussian beam), 5: X/Y factorizable beam (e.g., a rectangular LED emitter)
model.MC.beam.xFocus              = 0; % [cm] x position of focus
model.MC.beam.yFocus              = 0; % [cm] y position of focus
model.MC.beam.zFocus              = 0; % [cm] z position of focus
model.MC.beam.theta               = 0; % [rad] Polar angle of beam center axis
model.MC.beam.phi                 = 0; % [rad] Azimuthal angle of beam center axis

model.MC.useLightCollector        = true;
model.MC.LC.x                     = 0; % [cm] x position of either the center of the objective lens focal plane or the fiber tip
model.MC.LC.y                     = -0.05; % [cm] y position
model.MC.LC.z                     = 0.15; % [cm] z position

model.MC.LC.theta                 = 3*pi/4; % [rad] Polar angle of direction the light collector is facing
model.MC.LC.phi                   = pi/2; % [rad] Azimuthal angle of direction the light collector is facing

model.MC.LC.f                     = Inf; % [cm] Focal length of the objective lens (if light collector is a fiber, set this to Inf).
model.MC.LC.diam                  = .1; % [cm] Diameter of the light collector aperture. For an ideal thin lens, this is 2*f*tan(asin(NA)).
model.MC.LC.fieldSize             = .1; % [cm] Field Size of the imaging system (diameter of area in object plane that gets imaged). Only used for finite f.
model.MC.LC.NA                    = 0.22; % [-] Fiber NA. Only used for infinite f.

model.MC.LC.res                   = 1; % X and Y resolution of light collector in pixels, only used for finite f

%% Looping over the different tissue thicknesses
t_vec = linspace(0,0.1,21); % Thicknesses to simulate
power_vec = zeros(1,length(t_vec));
fprintf('%2d/%2d\n',0,length(t_vec));
for i=1:length(t_vec)
    fprintf('\b\b\b\b\b\b%2d/%2d\n',i,length(t_vec)); % Simple progress indicator
    
    % Adjust geometry
    model.G.geomFuncParams      = {t_vec(i)}; % Cell array containing any additional parameters to pass into the geometry function, such as media depths, inhomogeneity positions, radii etc.
    plot(model,'G');
    
    % Run MC
    model = runMonteCarlo(model);
    
    % Post-processing
    power_vec(i) = model.MC.LC.image; % "image" is in this case just a scalar, the normalized power collected by the fiber.
end
plot(model,'MC');

%% Plotting the collected power vs. tissue thickness
figure;clf;
plot(t_vec,power_vec,'Linewidth',2);
set(gcf,'Position',[40 80 800 550]);
xlabel('Slab thickness [cm]');
ylabel('Normalized power collected by fiber');
set(gca,'FontSize',18);grid on; grid minor;

%% Geometry function(s)
% A geometry function takes as input X,Y,Z matrices as returned by the
% "ndgrid" MATLAB function as well as any parameters the user may have
% provided in the definition of Ginput. It returns the media matrix M,
% containing numerical values indicating the media type (as defined in
% mediaPropertiesFunc) at each voxel location.
function M = geometryDefinition_VariableThicknessStandardTissue(X,Y,Z,parameters)
    M = ones(size(X)); % Air
    M(Z > 0.1 - parameters{1} & Z < 0.1) = 2; % "Standard" tissue
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