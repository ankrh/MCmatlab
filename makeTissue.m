function makeTissue
%% Define parameters (user-specified)
wavelength  = 850;      % [nm] set the range of wavelengths of the Monte Carlo simulation
name = 'hair';          % name of the simulation
nx = 100;               % number of bins in the x direction
ny = 100;                % number of bins in the y direction
nz = 50;                % number of bins in the z direction
Lx = 0.05;               % [cm] x size of simulation area
Ly = 0.05;               % [cm] y size of simulation area
Lz = 0.15;               % [cm] z size of simulation area

%% Calculate x,y,z vectors and grids
dx = Lx/nx;             % [cm] size of x bins
dy = Ly/ny;             % [cm] size of y bins
dz = Lz/nz;             % [cm] size of z bins
x  = ((0:nx-1)-(nx-1)/2)*dx;
y  = ((0:ny-1)-(ny-1)/2)*dy;
z  = ((0:nz-1)+1/2)*dz;
[X,Y,Z] = ndgrid(x,y,z);

%% Define tissue T(x,y,z) (user-specified)
% Dentin example:
% T = 3*ones(nx,ny,nz,'uint8'); % fill background with air
% T(Y.^2 + (Z-nz*dz/2).^2 < 0.300^2) = 2; % enamel
% T(Y.^2 + (Z-nz*dz/2).^2 < 0.170^2) = 1; % dentin
% T(Y.^2 + (Z-nz*dz/2).^2 < 0.050^2) = 5; % blood

% Hair example:
zsurf = 0.02;  % position of gel/skin surface[cm]
epd_thick = 0.01; % thickness of the epidermis [cm]
hair_radius = 0.0075/2; % diameter varies from 17 - 180 micrometers, should increase with colouring and age
hair_bulb_semiminor = 1.7*hair_radius; % [cm]
hair_bulb_semimajor = sqrt(2)*hair_bulb_semiminor;
hair_depth = 0.1; % varies from 0.06-0.3cm
papilla_semiminor = hair_bulb_semiminor*5/12;
papilla_semimajor = sqrt(2)*papilla_semiminor;

T = 4*ones(nx,ny,nz,'uint8'); % water (gel)
T(Z > zsurf) = 10; % epidermis
T(Z > zsurf+epd_thick) = 9; % dermis
T(X.^2 + Y.^2 < hair_radius^2 & Z < zsurf+hair_depth) = 15; % hair
T((X/hair_bulb_semiminor).^2 + (Y/hair_bulb_semiminor).^2 + ((Z-(zsurf+hair_depth))/hair_bulb_semimajor).^2 < 1) = 15; % hair
T((X/papilla_semiminor).^2 + (Y/papilla_semiminor).^2 + ((Z-(zsurf+hair_depth+hair_bulb_semimajor-papilla_semimajor))/papilla_semimajor).^2 < 1) = 9; % dermis (papilla)

%% Discard the unused tissue types
tissueList = makeTissueList(wavelength);
nT = length(unique(T)); % Number of different tissues in simulation
tissueMap = zeros(1,length(tissueList),'uint8');
tissueMap(unique(T)) = 1:nT;
tissueList = tissueList(unique(T));
T = tissueMap(T); % Reduced tissue matrix, using only numbers from 1 up to the number of tissues actually used

%% Make plots
figure(1);clf;
plotVolumetric(x,y,z,T,tissueList)
title('Tissue type illustration');
drawnow;

figure(7);clf;
plotTissueProperties(tissueList);
drawnow;

%% Save output
save(['.\Data\' name '.mat'],'x','y','z','tissueList','T','wavelength')
fprintf('.\\Data\\%s.mat saved\n',name)

