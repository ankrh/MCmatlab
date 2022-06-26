%% Description
% This example is another illustration of MC simulations inside a for loop,
% this time simulating a pencil beam incident on a 100µm slab of scattering
% medium with a variable (parametrically sweeped) scattering anisotropy g.
% g is passed in through the mediaPropParams field and used within
% mediaPropertiesFunc. Light is collected in transmission at a 45° angle in
% a fiber, similar to example 8. At the end of the script, collected power
% as a function of g is plotted. The power is seen to be zero for g = +- 1,
% which is because then the light can only be scattered exactly forward or
% backward. The max is at about 0.6, fitting well with a single scattering
% event at 45°. There is a secondary hump around -0.7, which fits with
% photons experiencing two scattering events at a scattering angle of
% 157.5°.
% 
% As in example 9, calcNFR is again set to false to speed up the simulation
% slightly.
% 
% The results of the last run will be plotted and MCmatlab will give a
% warning because in that last one, no photons were collected in the
% detector. This is not a problem in our case since we were expecting this.

%% Common MCmatlab abbreviations:
% G: Geometry, MC: Monte Carlo, FMC: Fluorescence Monte Carlo, HS: Heat
% simulation, M: Media array, LS: Light source, LC: Light collector, FPID:
% Focal plane intensity distribution, AID: Angular intensity distribution,
% NI: Normalized irradiance, NFR: Normalized fluence rate, FR: Fluence
% rate, FD: Fractional damage.

%% Geometry definition
model = MCmatlab.model;

model.G.silentMode        = true; % Disables command window text and progress indication

model.G.nx                = 21; % Number of bins in the x direction
model.G.ny                = 21; % Number of bins in the y direction
model.G.nz                = 21; % Number of bins in the z direction
model.G.Lx                = .1; % [cm] x size of simulation cuboid
model.G.Ly                = .1; % [cm] y size of simulation cuboid
model.G.Lz                = .01; % [cm] z size of simulation cuboid

model.G.mediaPropertiesFunc = @mediaPropertiesFunc; % Media properties defined as a function at the end of this file
model.G.geomFunc          = @geometryDefinition; % Function to use for defining the distribution of media in the cuboid. Defined at the end of this m file.

%% Monte Carlo simulation
model.MC.silentMode               = true; % Disables command window text and progress indication
model.MC.useAllCPUs               = true; % If false, MCmatlab will leave one processor unused. Useful for doing other work on the PC while simulations are running.
model.MC.simulationTimeRequested  = 2/60; % [min] Time duration of the simulation
model.MC.calcNFR                  = false; % (Default: true) If true, the 3D fluence rate output matrix NFR will be calculated. Set to false if you have a light collector and you're only interested in the image output.

model.MC.matchedInterfaces        = true; % Assumes all refractive indices are the same
model.MC.boundaryType             = 1; % 0: No escaping boundaries, 1: All cuboid boundaries are escaping, 2: Top cuboid boundary only is escaping, 3: Top and bottom boundaries are escaping, while the side boundaries are cyclic
model.MC.wavelength               = 532; % [nm] Excitation wavelength, used for determination of optical properties for excitation light

model.MC.LS.sourceType   = 0; % 0: Pencil beam, 1: Isotropically emitting line or point source, 2: Infinite plane wave, 3: Laguerre-Gaussian LG01 beam, 4: Radial-factorizable beam (e.g., a Gaussian beam), 5: X/Y factorizable beam (e.g., a rectangular LED emitter)
model.MC.LS.xFocus       = 0; % [cm] x position of focus
model.MC.LS.yFocus       = 0; % [cm] y position of focus
model.MC.LS.zFocus       = 0; % [cm] z position of focus
model.MC.LS.theta        = 0; % [rad] Polar angle of beam center axis
model.MC.LS.phi          = 0; % [rad] Azimuthal angle of beam center axis

model.MC.useLightCollector        = true;
model.MC.LC.x                     = 0; % [cm] x position of either the center of the objective lens focal plane or the fiber tip
model.MC.LC.y                     = -0.05; % [cm] y position
model.MC.LC.z                     = 0.06; % [cm] z position

model.MC.LC.theta                 = 3*pi/4; % [rad] Polar angle of direction the light collector is facing
model.MC.LC.phi                   = pi/2; % [rad] Azimuthal angle of direction the light collector is facing

model.MC.LC.f                     = Inf; % [cm] Focal length of the objective lens (if light collector is a fiber, set this to Inf).
model.MC.LC.diam                  = .1; % [cm] Diameter of the light collector aperture. For an ideal thin lens, this is 2*f*tan(asin(NA)).
model.MC.LC.fieldSize             = .1; % [cm] Field Size of the imaging system (diameter of area in object plane that gets imaged). Only used for finite f.
model.MC.LC.NA                    = 0.22; % [-] Fiber NA. Only used for infinite f.

model.MC.LC.res                   = 1; % X and Y resolution of light collector in pixels, only used for finite f

%% Looping over the different scattering anisotropies g
g_vec = linspace(-1,1,21); % g values to simulate
power_vec = zeros(1,length(g_vec));
fprintf('%2d/%2d\n',0,length(g_vec));
for i=1:length(g_vec)
    fprintf('\b\b\b\b\b\b%2d/%2d\n',i,length(g_vec)); % Simple progress indicator
    
    % Adjust media properties
    model.G.mediaPropParams   = {g_vec(i)}; % Cell array containing any additional parameters to be passed to the mediaPropertiesFunc function
    
    % Run MC
    model = runMonteCarlo(model);
    
    % Post-processing
    power_vec(i) = model.MC.LC.image; % "image" is in this case just a scalar, the normalized power collected by the fiber.
end
plot(model,'G');
plot(model,'MC');

%% Plotting the collected power vs. scattering anisotropy g
figure;clf;
plot(g_vec,power_vec,'Linewidth',2);
set(gcf,'Position',[40 80 800 550]);
xlabel('Scattering anisotropy g');
ylabel('Normalized power collected by fiber');
set(gca,'FontSize',18);grid on; grid minor;

%% Geometry function(s)
% A geometry function takes as input X,Y,Z matrices as returned by the
% "ndgrid" MATLAB function as well as any parameters the user may have
% provided in the definition of Ginput. It returns the media matrix M,
% containing numerical values indicating the media type (as defined in
% mediaPropertiesFunc) at each voxel location.
function M = geometryDefinition(X,Y,Z,parameters)
    M = ones(size(X)); % Variable g medium
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
    mediaProperties(j).name  = 'variable g medium';
    mediaProperties(j).mua   = 10; % [cm^-1]
    mediaProperties(j).mus   = 100; % [cm^-1]
    mediaProperties(j).g = parameters{1};
end
