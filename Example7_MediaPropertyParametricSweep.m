addpath([fileparts(mfilename('fullpath')) '/helperfuncs']); % The helperfuncs folder is added to the path for the duration of this MATLAB session

%% Description
% This example is another illustration of MC simulations inside a for loop,
% this time simulating a pencil beam incident on a 100µm slab of scattering
% medium with a variable (parametrically sweeped) scattering anisotropy g.
% g is passed in through the mediaPropParams field and used within
% getMediaProperties. Light is collected in transmission at a 45° angle in
% a fiber, similar to example 6. At the end of the script, collected power
% as a function of g is plotted. The power is seen to be zero for g = +- 1,
% which is because then the light can only be scattered exactly forward or
% backward. The max is at about 0.6, fitting well with a single scattering
% event at 45°. There is a secondary hump around -0.7, which fits with
% photons experiencing two scattering events at a scattering angle of
% 157.5°.

g_vec = linspace(-1,1,21); % g values to simulate
power_vec = zeros(1,length(g_vec));
fprintf('%2d/%2d\n',0,length(g_vec));
for i=1:length(g_vec)
fprintf('\b\b\b\b\b\b%2d/%2d\n',i,length(g_vec)); % Simple progress indicator
%% Geometry definition
clear Ginput
Ginput.silentMode        = true; % Disables command window text and progress indication
Ginput.matchedInterfaces = true; % Assumes all refractive indices are 1
Ginput.boundaryType      = 1; % 0: No boundaries, 1: All cuboid boundaries, 2: Top cuboid boundary only

Ginput.wavelength        = 532; % [nm] Excitation wavelength, used for determination of optical properties for excitation light
Ginput.mediaPropParams   = {g_vec(i)}; % Cell array containing any additional parameters to be passed to the getMediaProperties function

Ginput.nx                = 21; % Number of bins in the x direction
Ginput.ny                = 21; % Number of bins in the y direction
Ginput.nz                = 21; % Number of bins in the z direction
Ginput.Lx                = .1; % [cm] x size of simulation cuboid
Ginput.Ly                = .1; % [cm] y size of simulation cuboid
Ginput.Lz                = .01; % [cm] z size of simulation cuboid

Ginput.GeomFunc          = @GeometryDefinition_MediaPropertyParametricSweep; % Function to use for defining the distribution of media in the cuboid. Defined at the end of this m file.

% Execution, do not modify the next two lines:
Goutput = defineGeometry(Ginput);
% plotMCmatlabGeom(Goutput);

%% Monte Carlo simulation
clear MCinput
MCinput.silentMode               = true; % Disables command window text and progress indication
MCinput.useAllCPUs               = true; % If false, MCmatlab will leave one processor unused. Useful for doing other work on the PC while simulations are running.
MCinput.simulationTime           = 2/60; % [min] Time duration of the simulation

MCinput.Beam.beamType            = 0; % 0: Pencil beam, 1: Isotropically emitting point source, 2: Infinite plane wave, 3: Gaussian focus, Gaussian far field beam, 4: Gaussian focus, top-hat far field beam, 5: Top-hat focus, Gaussian far field beam, 6: Top-hat focus, top-hat far field beam, 7: Laguerre-Gaussian LG01 beam
MCinput.Beam.xFocus              = 0; % [cm] x position of focus
MCinput.Beam.yFocus              = 0; % [cm] y position of focus
MCinput.Beam.zFocus              = 0; % [cm] z position of focus
MCinput.Beam.theta               = 0; % [rad] Polar angle of beam center axis
MCinput.Beam.phi                 = 0; % [rad] Azimuthal angle of beam center axis
MCinput.Beam.waist               = 0.005; % [cm] Beam waist 1/e^2 radius
MCinput.Beam.divergence          = 5/180*pi; % [rad] Beam divergence 1/e^2 half-angle of beam (for a diffraction limited Gaussian beam, this is G.wavelength*1e-9/(pi*MCinput.Beam.waist*1e-2))

MCinput.LightCollector.x         = 0; % [cm] x position of either the center of the objective lens focal plane or the fiber tip
MCinput.LightCollector.y         = -0.05; % [cm] y position
MCinput.LightCollector.z         = 0.06; % [cm] z position

MCinput.LightCollector.theta     = 3*pi/4; % [rad] Polar angle of direction the light collector is facing
MCinput.LightCollector.phi       = pi/2; % [rad] Azimuthal angle of direction the light collector is facing

MCinput.LightCollector.f         = Inf; % [cm] Focal length of the objective lens (if light collector is a fiber, set this to Inf).
MCinput.LightCollector.diam      = .1; % [cm] Diameter of the light collector aperture. For an ideal thin lens, this is 2*f*tan(asin(NA)).
MCinput.LightCollector.FieldSize = .1; % [cm] Field Size of the imaging system (diameter of area in object plane that gets imaged). Only used for finite f.
MCinput.LightCollector.NA        = 0.22; % [-] Fiber NA. Only used for infinite f.

MCinput.LightCollector.res       = 1; % X and Y resolution of light collector in pixels, only used for finite f

% Execution, do not modify the next three lines:
MCinput.G = Goutput;
MCoutput = runMonteCarlo(MCinput);
% plotMCmatlab(MCinput,MCoutput);

%% Post-processing
power_vec(i) = MCoutput.Image; % "Image" is in this case just a scalar, the normalized power collected by the fiber.
end

figure;clf;
plot(g_vec,power_vec,'Linewidth',2);
set(gcf,'Position',[40 80 1100 650]);
xlabel('Scattering anisotropy g');
ylabel('Normalized power collected by fiber');
set(gca,'FontSize',18);grid on; grid minor;

%% Geometry function(s)
% A geometry function takes as input X,Y,Z matrices as returned by the
% "ndgrid" MATLAB function as well as any parameters the user may have
% provided in the definition of Ginput. It returns the media matrix M,
% containing numerical values indicating the media type (as defined in
% getMediaProperties) at each voxel location.
function M = GeometryDefinition_MediaPropertyParametricSweep(X,Y,Z,parameters)
M = 21*ones(size(X)); % Variable g medium
end
