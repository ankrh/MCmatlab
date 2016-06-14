% Generates plots of the used tissue properties
clear all
close all 
clc
%% Load spectral library
load SpectralLIB    
LOADED = 1;
nm = 300:5:1000;
%% Create tissueList
% dermis
B = 0.002; 
S = 0.67;
W = 0.65;
M = 0;
musp500 = 42.4;
fray    = 0.62;
bmie    = 1.0;
X = [B*S B*(1-S) W M]';
for j=1:length(nm)
MU(:,1) = interp1(nmLIB,muaoxy,nm(j));
MU(:,2) = interp1(nmLIB,muadeoxy,nm(j));
MU(:,3) = interp1(nmLIB,muawater,nm(j));
MU(:,4) = interp1(nmLIB,muamel,nm(j));
mua_dermis(j) = MU*X;
end
musp_dermis = musp500*(fray*(nm/500).^-4 + (1-fray)*(nm/500).^-bmie);
gg = 0.90;
mus = musp_dermis/(1-gg);
g   = gg;
HC =3391*1109e-6;
%D(j) =1109e-6;
TC =0.37e-2;

% epidermis
B = 0;
S = 0.75;
W = 0.75;
M = 0.03;
musp500 = 40;
fray    = 0.0;
bmie    = 1.0;
X = [B*S B*(1-S) W M]';
for j=1:length(nm)
    MU(:,1) = interp1(nmLIB,muaoxy,nm(j));
    MU(:,2) = interp1(nmLIB,muadeoxy,nm(j));
    MU(:,3) = interp1(nmLIB,muawater,nm(j));
    MU(:,4) = interp1(nmLIB,muamel,nm(j));
    mua_epidermis(j) = MU*X;
end
musp_epidermis = musp500*(fray*(nm/500).^-4 + (1-fray)*(nm/500).^-bmie);
gg = 0.90;
mus = musp_epidermis/(1-gg);
g   = gg;
HC =3391*1109e-6;
%D(j) =1109e-6;
TC =0.37e-2;

% Blood
B       = 1.00;
S       = 0.75;
W       = 0.95;
M       = 0;
musp500 = 10;
fray    = 0.0;
bmie    = 1.0;
gg      = 0.90;
X = [B*S B*(1-S) W M]';
for j=1:length(nm)
    MU(:,1) = interp1(nmLIB,muaoxy,nm(j));
    MU(:,2) = interp1(nmLIB,muadeoxy,nm(j));
    MU(:,3) = interp1(nmLIB,muawater,nm(j));
    MU(:,4) = interp1(nmLIB,muamel,nm(j));
    mua_blood(j) = MU*X;
end
musp_blood = musp500*(fray*(nm/500).^-4 + (1-fray)*(nm/500).^-bmie);
mus = musp_blood/(1-gg);
g   = gg;
HC =3617*1050e-6;
%D(j) =1050e-6;
TC =0.52e-2;

% Hair
musp_hair = 9.6e10./nm.^3; % from Bashkatov 2002, for a dark hair

figure
semilogy(nm,mua_dermis,'-k',nm,mua_epidermis,'--b',nm,mua_blood,'-.r','Linewidth',3)
set(gca,'FontSize',24);xlabel('Wavelength [nm]');ylabel('Absoption Coefficient')
legend('Dermis','Epidermis','Blood')

figure
semilogy(nm,musp_dermis,'-k',nm,musp_epidermis,'--b',nm,musp_blood,'-.r',nm,musp_hair,':m','Linewidth',3)
set(gca,'FontSize',24);xlabel('Wavelength [nm]');ylabel('Reduced Scattering Coefficient')
legend('Dermis','Epidermis','Blood','Hair')


semilogy(nm,musp_hair)