function tissueList = makeTissueList(nm)
%function tissueProps = makeTissueList(nm)
%   Returns the tissue optical properties at the wavelength nm: and the
%   thermal properties, note that the data for skull and brainmatter are
%   not correct, they are simply placeholders
%       tissueProps = [mua; mus; g; HC; TC]'; HC is heat capacity
%       [J/cm^3/K], D is density kg/cm^3, TC is thermal conductivity W/cm/K
%   Uses 
%       SpectralLIB.mat

%% Load spectral library
load spectralLIB.mat
%   muadeoxy      701x1              5608  double              
%   muamel        701x1              5608  double              
%   muaoxy        701x1              5608  double              
%   muawater      701x1              5608  double              
%   musp          701x1              5608  double              
%   nmLIB         701x1              5608  double              
MU(:,1) = interp1(nmLIB,muaoxy,nm);
MU(:,2) = interp1(nmLIB,muadeoxy,nm);
MU(:,3) = interp1(nmLIB,muawater,nm);
MU(:,4) = interp1(nmLIB,muamel,nm);

%% Create tissueList

j=1;
tissueList(j).name  = 'escape';
tissueList(j).mua   = 0.0001;
tissueList(j).mus   = 1.0;
tissueList(j).g     = 1.0;
tissueList(j).HC    = 1003*1e-6; % Heat Capacity
tissueList(j).D     = 1e-6; % Density, [kg/cm^3]
tissueList(j).TC    = 3e-4; % Thermal Conductivity

j=2;
tissueList(j).name  = 'air';
tissueList(j).mua   = 0.001;
tissueList(j).mus   = 10;
tissueList(j).g     = 1.0;
tissueList(j).HC    = 1003*1e-6;
tissueList(j).D     = 1e-6;
tissueList(j).TC    = 3e-4;

j=3;
tissueList(j).name  = 'blood';
B       = 1.00;
S       = 0.75;
W       = 0.95;
M       = 0;
musp500 = 10;
fray    = 0.0;
bmie    = 1.0;
gg      = 0.90;
musp = musp500*(fray*(nm/500).^-4 + (1-fray)*(nm/500).^-bmie);
X = [B*S B*(1-S) W M]';
tissueList(j).mua = MU*X;
tissueList(j).mus = musp/(1-gg);
tissueList(j).g   = gg;
tissueList(j).HC    = 3617*1050e-6;
tissueList(j).D     = 1050e-6;
tissueList(j).TC    = 0.52e-2;

j = 4;
tissueList(j).name = 'dermis';
B = 0.002; 
S = 0.67;
W = 0.65;
M = 0;
musp500 = 42.4;
fray    = 0.62;
bmie    = 1.0;
gg      = 0.90;
musp = musp500*(fray*(nm/500).^-4 + (1-fray)*(nm/500).^-bmie);
X = [B*S B*(1-S) W M]';
tissueList(j).mua = MU*X;
tissueList(j).mus = musp/(1-gg);
tissueList(j).g   = gg;
tissueList(j).HC    = 3391*1109e-6;
tissueList(j).D     = 1109e-6;
tissueList(j).TC    = 0.37e-2;

j=5;
tissueList(j).name  = 'epidermis';
B = 0;
S = 0.75;
W = 0.75;
M = 0.03;
musp500 = 40;
fray    = 0.0;
bmie    = 1.0;
gg      = 0.90;
musp = musp500*(fray*(nm/500).^-4 + (1-fray)*(nm/500).^-bmie);
X = [B*S B*(1-S) W M]';
tissueList(j).mua = MU*X;
tissueList(j).mus = musp/(1-gg);
tissueList(j).g   = gg;
tissueList(j).HC    = 3391*1109e-6;
tissueList(j).D     = 1109e-6;
tissueList(j).TC    = 0.37e-2;

j=6;
tissueList(j).name  = 'skull';
B = 0.0005;
S = 0.75;
W = 0.35;
M = 0;
musp500 = 30;
fray    = 0.0;
bmie    = 1.0;
gg      = 0.90;
musp = musp500*(fray*(nm/500).^-4 + (1-fray)*(nm/500).^-bmie);
X = [B*S B*(1-S) W M]';
tissueList(j).mua = MU*X;
tissueList(j).mus = musp/(1-gg);
tissueList(j).g   = gg;
tissueList(j).HC    = 3391*1109e-6;
tissueList(j).D     = 1109e-6;
tissueList(j).TC    = 0.37e-2;

j=7;
tissueList(j).name = 'gray matter';
B = 0.01;
S = 0.75;
W = 0.75;
M = 0;
musp500 = 20;
fray    = 0.2;
bmie    = 1.0;
gg      = 0.90;
musp = musp500*(fray*(nm/500).^-4 + (1-fray)*(nm/500).^-bmie);
X = [B*S B*(1-S) W M]';
tissueList(j).mua = MU*X;
tissueList(j).mus = musp/(1-gg);
tissueList(j).g   = gg;
tissueList(j).HC    = 3391*1109e-6;
tissueList(j).D     = 1109e-6;
tissueList(j).TC    = 0.37e-2;

j=8;
tissueList(j).name  = 'white matter';
B = 0.01;
S = 0.75;
W = 0.75;
M = 0;
musp500 = 20;
fray    = 0.2;
bmie    = 1.0;
gg      = 0.90;
musp = musp500*(fray*(nm/500).^-4 + (1-fray)*(nm/500).^-bmie);
X = [B*S B*(1-S) W M]';
tissueList(j).mua = MU*X;
tissueList(j).mus = musp/(1-gg);
tissueList(j).g   = gg;
tissueList(j).HC    = 3391*1109e-6;
tissueList(j).D     = 1109e-6;
tissueList(j).TC    = 0.37e-2;

j=9;
tissueList(j).name = 'Hair';
B = 0;
S = 0.75;
W = 0.75;
M = 0.03;
X = 2*[B*S B*(1-S) W M]'; 
musp_hair = 9.6e10./nm.^3; % from Bashkatov 2002, for a dark hair
tissueList(j).mua   = MU*X; %Here the hair is set to absorb twice as much as the epidermis
tissueList(j).mus   = musp_hair/(1-0.9); %% use (1-g) for the denominator
tissueList(j).g     = 0.9;
tissueList(j).HC    = 1.53*1.3; % Thermal data has been approximated using the data for horn, as horn and hair are both composed of keratin
tissueList(j).D     = 1.3e-3;
tissueList(j).TC    = 6.3e-3;

j=10;
tissueList(j).name = 'Water'; % Used for the gel
% mus500, fray and bmie for water have been found by fitting to the data
% from Buiteveld, H., JHM Hakvoort, and M. Donze. 
% THE OPTICAL PROPERTIES OF PURE WATER. OCEAN OPTICS XII 2258 (1994):
% 174-183.
mus500 = 0.2102; %[cm^{-1}]
fray    = 1.122;
bmie    = 2.482;
tissueList(j).mua   = MU(:,3);
tissueList(j).mus   = 10+mus500*(fray*(nm/500).^-4 + (1-fray)*(nm/500).^-bmie);
tissueList(j).g     = 0.99;
tissueList(j).HC    = 4.1796; %http://en.wikipedia.org/wiki/Heat_capacity
tissueList(j).D     = 1e-3;
tissueList(j).TC    = 0.6e-2; %http://en.wikipedia.org/wiki/List_of_thermal_conductivities 14-01-2014

j=11;
tissueList(j).name  = 'standard tissue';
tissueList(j).mua   = 1;
tissueList(j).mus   = 100;
tissueList(j).g     = 0.90;
tissueList(j).HC    = 4.1796; %http://en.wikipedia.org/wiki/Heat_capacity
tissueList(j).D     = 1e-3;
tissueList(j).TC    = 0.6e-2; %http://en.wikipedia.org/wiki/List_of_thermal_conductivities 14-01-2014

fprintf('---- tissueList ------ \tmua   \tmus  \tg  \tmusp\n')
for i=1:length(tissueList)
    fprintf('%d\t%15s\t%0.4f\t%0.1f\t%0.3f\t%0.1f\n\n',...
        i,tissueList(i).name, tissueList(i).mua,tissueList(i).mus,tissueList(i).g,...
        tissueList(i).mus*(1-tissueList(i).g))
end

