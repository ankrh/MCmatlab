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
tissueList(j).name  = 'dentin';
tissueList(j).mua   = 0.0001;
tissueList(j).mus   = 1.0;
tissueList(j).g     = 0.93;
tissueList(j).HC    = 1.59; % Heat Capacity
tissueList(j).D     = 1e-6; % Density, [kg/cm^3]
tissueList(j).TC    = 5.69e-3; % Thermal Conductivity


j=2;
tissueList(j).name  = 'enamel';
tissueList(j).mua   = 0.001;
tissueList(j).mus   = 10;
tissueList(j).g     = 0.96;
tissueList(j).HC    = 0.71;
tissueList(j).D     = 1e-6;
tissueList(j).TC    = 9.34e-3;

j=3;
tissueList(j).name  = 'air';
tissueList(j).mua   = 0.001;
tissueList(j).mus   = 10;
tissueList(j).g     = 1.0;
tissueList(j).HC    = 1003*1e-6;
tissueList(j).D     = 1;
tissueList(j).TC    = 3e-4;

j=4;
tissueList(j).name  = 'water';
tissueList(j).mua   = 0.001;
tissueList(j).mus   = 10;
tissueList(j).g     = 1.0;
tissueList(j).HC    = 4.2;
tissueList(j).D     = 1.0;
tissueList(j).TC    = 6.2; % between 0.561 and 0.679 W/m/K in cm its 6.2

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
tissueList(j).HC    = 3617*1050e-6;
tissueList(j).D     = 1050e-6;
tissueList(j).TC    = 0.52e-2;




fprintf('---- tissueList ------ \tmua   \tmus  \tg  \tmusp\n')
for i=1:length(tissueList)
    fprintf('%d\t%15s\t%0.4f\t%0.1f\t%0.3f\t%0.1f\n\n',...
        i,tissueList(i).name, tissueList(i).mua,tissueList(i).mus,tissueList(i).g,...
        tissueList(i).mus*(1-tissueList(i).g))
end

