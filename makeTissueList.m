function tissueProps = makeTissueList(nm)
%function tissueProps = makeTissueList(nm)
%   Returns the tissue optical properties at the wavelength nm: and the
%   thermal properties, note that the data for skull and brainmatter are
%   not correct, they are simply placeholders
%       tissueProps = [mua; mus; g; HC; TC]'; HC is heat capacity
%       [J/cm^3/K], D is density kg/cm^3, TC is thermal conductivity W/cm/K
%   Uses 
%       SpectralLIB.mat
global tissue

%% Load spectral library
load SpectralLIB
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
LOADED = 1;

%% Create tissueList

j=1;
tissue(j).s = 'escape';
mua(j) = 0.0001;
mus(j) = 1.0;
g(j)   = 1.0;
HC(j) =1003*1e-6;
%D(j) =1e-6;
TC(j) =3e-4;

j=2;
tissue(j).s = 'air';
mua(j) = 0.1;
mus(j) = 10;
g(j)   = 1.0;
HC(j) =1003*1e-6;
%D(j) =1e-6;
TC(j) =3e-4;

j=3;
tissue(j).s = 'blood';
B       = 1.00;
S       = 0.75;
W       = 0.95;
M       = 0;
musp500 = 10;
fray    = 0.0;
bmie    = 1.0;
gg      = 0.90;
X = [B*S B*(1-S) W M]';
mua(j) = MU*X;
musp = musp500*(fray*(nm/500).^-4 + (1-fray)*(nm/500).^-bmie);
mus(j) = musp/(1-gg);
g(j)   = gg;
HC(j) =3617*1050e-6;
%D(j) =1050e-6;
TC(j) =0.52e-2;

j = 4;
tissue(j).s = 'dermis';
B = 0.002; 
S = 0.67;
W = 0.65;
M = 0;
musp500 = 42.4;
fray    = 0.62;
bmie    = 1.0;
X = [B*S B*(1-S) W M]';
mua(j) = MU*X;
musp = musp500*(fray*(nm/500).^-4 + (1-fray)*(nm/500).^-bmie);
gg = 0.90;
mus(j) = musp/(1-gg);
g(j)   = gg;
HC(j) =3391*1109e-6;
%D(j) =1109e-6;
TC(j) =0.37e-2;

j=5;
tissue(j).s = 'epidermis';
B = 0;
S = 0.75;
W = 0.75;
M = 0.03;
musp500 = 40;
fray    = 0.0;
bmie    = 1.0;
X = [B*S B*(1-S) W M]';
mua(j) = MU*X;
musp = musp500*(fray*(nm/500).^-4 + (1-fray)*(nm/500).^-bmie);
gg = 0.90;
mus(j) = musp/(1-gg);
g(j)   = gg;
HC(j) =3391*1109e-6;
%D(j) =1109e-6;
TC(j) =0.37e-2;

j=6;
tissue(j).s = 'skull';
B = 0.0005;
S = 0.75;
W = 0.35;
M = 0;
musp500 = 30;
fray    = 0.0;
bmie    = 1.0;
X = [B*S B*(1-S) W M]';
mua(j) = MU*X;
musp = musp500*(fray*(nm/500).^-4 + (1-fray)*(nm/500).^-bmie);
gg = 0.90;
mus(j) = musp/(1-gg);
g(j)   = gg;
HC(j) =3391*1109e-6;
%D(j) =1109e-6;
TC(j) =0.37e-2;

j=7;
tissue(j).s = 'gray matter';
B = 0.01;
S = 0.75;
W = 0.75;
M = 0;
musp500 = 20;
fray    = 0.2;
bmie    = 1.0;
X = [B*S B*(1-S) W M]';
mua(j) = MU*X;
musp = musp500*(fray*(nm/500).^-4 + (1-fray)*(nm/500).^-bmie);
gg = 0.90;
mus(j) = musp/(1-gg);
g(j)   = gg;
HC(j) =3391*1109e-6;
%D(j) =1109e-6;
TC(j) =0.37e-2;

j=8;
tissue(j).s = 'white matter';
B = 0.01;
S = 0.75;
W = 0.75;
M = 0;
musp500 = 20;
fray    = 0.2;
bmie    = 1.0;
X = [B*S B*(1-S) W M]';
mua(j) = MU*X;
musp = musp500*(fray*(nm/500).^-4 + (1-fray)*(nm/500).^-bmie);
gg = 0.90;
mus(j) = musp/(1-gg);
g(j)   = gg;
HC(j) =3391*1109e-6;
%D(j) =1109e-6;
TC(j) =0.37e-2;

j=9;
tissue(j).s = 'Hair';
B = 0;
S = 0.75;
W = 0.75;
M = 0.03;
X = 2*[B*S B*(1-S) W M]'; 
mua(j) = MU*X; %Here the hair is set to absorb twice as much as the epidermis
musp_hair = 9.6e10./nm.^3; % from Bashkatov 2002, for a dark hair
g(j)   = 0.9;
mus(j) = musp_hair/(1-g(j));
HC(j) =1.53*1.3; % Thermal data has been approximated using the data for horn, as horn and hair are both composed of keratin
%D(j) =1.3; [g/cm^{3}]
TC(j) =6.3e-3;

j=10;
tissue(j).s = 'Water'; % Used for the gel
% mus500, fray and bmie for water have been found by fitting to the data
% from Buiteveld, H., JHM Hakvoort, and M. Donze. 
% “THE OPTICAL PROPERTIES OF PURE WATER”. OCEAN OPTICS XII 2258 (1994):
% 174-183.
mus500 = 0.2102; %[cm^{-1}]
fray    = 1.122;
bmie    = 2.482;
mua(j)=MU(:,3);
mus(j)= 10+mus500*(fray*(nm/500).^-4 + (1-fray)*(nm/500).^-bmie);
g(j)=.99;
HC(j) =4.1796; %http://en.wikipedia.org/wiki/Heat_capacity
%D(j) =1e-6;
TC(j) =0.6e-2; %http://en.wikipedia.org/wiki/List_of_thermal_conductivities 14-01-2014

tissueProps = [mua; mus; g; HC; TC]';
disp(sprintf('---- tissueList ------ \tmua   \tmus  \tg  \tmusp'))
for i=1:length(mua)
    disp(sprintf('%d\t%15s\t\t%0.4f\t%0.2f\t%0.3f\t%0.1f',i,tissue(i).s, mua(i),mus(i),g(i),mus(i)*(1-g(i))))
end
disp(' ')

