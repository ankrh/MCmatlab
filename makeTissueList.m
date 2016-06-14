function tissueProps = makeTissueList(nm)
%function tissueProps = makeTissueList(nm)
%   Returns the tissue optical properties at the wavelength nm:
%       tissueProps = [mua; mus; g]';
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

j=2;
tissue(j).s = 'air';
mua(j) = 0.1;
mus(j) = 10;
g(j)   = 1.0;

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


tissueProps = [mua; mus; g]';
disp(sprintf('---- tissueList ------ \tmua   \tmus  \tg  \tmusp'))
for i=1:length(mua)
    disp(sprintf('%d\t%15s\t%0.4f\t%0.1f\t%0.3f\t%0.1f',i,tissue(i).s, mua(i),mus(i),g(i),mus(i)*(1-g(i))))
end
disp(' ')

