function tissueList = makeTissueList(nm)
%   Returns the tissue optical properties at the wavelength nm: and the thermal properties.
%   Note that the data for skull and brainmatter are not correct, they are simply placeholders
%
%   tissueList contains [mua; mus; g; VHC; D; TC] for each tissue type;
%       VHC is volumetric heat capacity [J/cm^3/K], D is density kg/cm^3, TC is thermal conductivity W/cm/K
%
%   Requires
%       SpectralLIB.mat

%% Updates
%   Anders K. Hansen & Dominik Marti, June 2017
%   Steven L. Jacques, August 2014

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
tissueList(j).name  = 'dentin';
tissueList(j).mua   = 0.04; %in cm ^ -1
tissueList(j).mus   = 270; %range between 260-280
tissueList(j).g     = 0.93;
tissueList(j).VHC   = 1260*2.2e-3; % Volumetric Heat Capacity [J/(cm^3*K)]
tissueList(j).D     = 2.200e-3; % Density, [kg/cm^3]
tissueList(j).TC    = 6e-3; % Thermal Conductivity [W/(cm*K)]

j=2;
tissueList(j).name = 'enamel';
tissueList(j).mua   = 0.01;
tissueList(j).mus   = 40;
tissueList(j).g     = 0.96;
tissueList(j).VHC   = 750*2.9e-3;
tissueList(j).D     = 2.900e-3;
tissueList(j).TC    = 9e-3;

j=3;
tissueList(j).name  = 'air';
tissueList(j).mua   = 1e-8;
tissueList(j).mus   = 1e-8;
tissueList(j).g     = 0;
tissueList(j).VHC   = 1.2e-3;
tissueList(j).D     = 1.2e-6;
tissueList(j).TC    = 0; % Real value is 2.6e-4, but we set it to zero to neglect the heat transport to air

j=4;
tissueList(j).name  = 'water';
tissueList(j).mua   = 0.001;
tissueList(j).mus   = 10;
tissueList(j).g     = 1.0;
tissueList(j).VHC   = 4.19;
tissueList(j).D     = 1e-3;
tissueList(j).TC    = 5.8e-3;

j=5;
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
tissueList(j).VHC   = 3617*1.050e-3;
tissueList(j).D     = 1.050e-3;
tissueList(j).TC    = 0.52e-2;

j=6;
tissueList(j).name  = 'fiber';
tissueList(j).mua   = 0.0001;
tissueList(j).mus   = 0.6666;
tissueList(j).g     = 0;
tissueList(j).VHC   = 703*2.203e-3;
tissueList(j).D     = 2.203e-3;
tissueList(j).TC    = 13.8e-3;

j=7;
tissueList(j).name  = 'vessel';
tissueList(j).mua   = 0.8;
tissueList(j).mus   = 230;
tissueList(j).g     = 0.9;
tissueList(j).VHC   = 4200*1.06e-3;
tissueList(j).D     = 1.06-3;
tissueList(j).TC    = 6.1e-3;

j=8;
tissueList(j).name  = 'patch';
tissueList(j).mua   = 1119;
tissueList(j).mus   = 15;
tissueList(j).g     = 0.8;
tissueList(j).VHC   = 5.363*1.048e-3;
tissueList(j).D     = 1.048e-3;
tissueList(j).TC    = 4.6e-3;

j=9;
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
tissueList(j).VHC =3391*1.109e-3;
tissueList(j).D =1.109e-3;
tissueList(j).TC =0.37e-2;

j=10;
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
tissueList(j).VHC =3391*1.109e-3;
tissueList(j).D =1.109e-3;
tissueList(j).TC =0.37e-2;

j=11;
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
tissueList(j).VHC = 3391*1.109e-3;
tissueList(j).D  = 1.109e-3;
tissueList(j).TC = 0.37e-2;

j=12;
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
tissueList(j).VHC = 3391*1.109e-3;
tissueList(j).D  = 1.109e-3;
tissueList(j).TC = 0.37e-2;

j=13;
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
tissueList(j).VHC = 3391*1.109e-3;
tissueList(j).D  = 1.109e-3;
tissueList(j).TC = 0.37e-2;

j=14;
tissueList(j).name  = 'standard tissue';
tissueList(j).mua   = 1;
tissueList(j).mus   = 100;
tissueList(j).g     = 0.90;
tissueList(j).VHC = 3391*1.109e-3;
tissueList(j).D  = 1.109e-3;
tissueList(j).TC = 0.37e-2;

j=15;
tissueList(j).name = 'hair';
B = 0;
S = 0.75;
W = 0.75;
M = 0.03;
gg = 0.9;
musp = 9.6e10./nm.^3; % from Bashkatov 2002, for a dark hair
X = 2*[B*S B*(1-S) W M]';
tissueList(j).mua = MU*X; %Here the hair is set to absorb twice as much as the epidermis
tissueList(j).mus = musp/(1-gg);
tissueList(j).g   = gg;
tissueList(j).VHC = 1530*1.3e-3; % Thermal data has been approximated using the data for horn, as horn and hair are both composed of keratin
tissueList(j).D   = 1.3e-3;
tissueList(j).TC  = 6.3e-3;
