function lookmcxyz

% lookmcxyz.m
%   Looks at myname_F.bin, created by mcxyz 
%   where myname is the name of the run: myname_T.bin, myname_H.mci
%
% Displays
%   Volumetric plot of tissue structure
%   Volumetric plot of fluence rate (intensity)
%   Volumetric plot of absorbed power per unit volume
%
% Saves
%   Fzy_data4.mat = Fzy y z zzs Fdet
%       Fzy(400,400,8) = 8 z,y images
%       Fdet(8,1) = signal [1/cm^2] @ detector fiber
%

%% USER CHOICES <---------- you must specify -----
directoryPath = 'exec/';
myname = 'dentin_sim_850';
wavelength = 850;
saveon_HeatSim = 1;

%% Load header file
H_mci = reportHmci(directoryPath,myname);

format compact

%% Load Fluence rate F(x,y,z) 
filename = sprintf('%s%s_F.bin',directoryPath,myname);
disp(['loading ' filename])
tic
    fid = fopen(filename, 'rb');
    data = fread(fid, H_mci.Nx*H_mci.Ny*H_mci.Nz, 'float');
    fclose(fid);
toc
F = reshape(data,H_mci.Nx,H_mci.Ny,H_mci.Nz); % F(x,y,z)

%% Load tissue structure in voxels, T(x,y,z) 
filename = sprintf('%s%s_T.bin',directoryPath,myname);
disp(['loading ' filename])
tic
    fid = fopen(filename, 'rb');
    data = fread(fid, H_mci.Nx*H_mci.Ny*H_mci.Nz, 'uint8=>uint8');
    fclose(fid);
toc
T = reshape(data,H_mci.Nx,H_mci.Ny,H_mci.Nz); % T(x,y,z)
clear data

%% Voxel center positions
x  = ((0:H_mci.Nx-1)-(H_mci.Nx-1)/2)*H_mci.dx;
y  = ((0:H_mci.Ny-1)-(H_mci.Ny-1)/2)*H_mci.dy;
z  = ((0:H_mci.Nz-1)+1/2)*H_mci.dz;
tissueList = makeTissueList(wavelength);

%% Make volumetric tissue plot

figure(1);clf;
plotVolumetric(x,y,z,T,tissueList);
title('Tissue type illustration');

%% Make volumetric fluence rate plot

figure(2);clf;
plotVolumetric(x,y,z,F);
title('Fluence rate (Intensity) [W/cm^2/W.delivered] ')

%% calculate power absorption

% F, matrix with fluency / input power
% T, Matrix containing tissue types
% tissueList, list of tissue properties tissueProps( tissue type, #) # = 1 is mua, # = 2 is mus, # = 3 is g

Ap = zeros(size(T));
for tissueNumber=1:length(tissueList)
   Ap(T==tissueNumber) = tissueList(tissueNumber).mua;
end

% test H_mci.dx*H_mci.dy*H_mci.dz*sum(sum(sum(Ap)))<=1, to make sure that less than 1W is
% absorbed per 1W of incident power, this seems to be the case

%% Make volumetric power absorption plot

figure(3);clf;
plotVolumetric(x,y,z,Ap.*F);
title('Absorbed power per unit volume [W/cm^3/W.delivered] ')

%% Save HeatSim input

if saveon_HeatSim==1
    clearvars F Azy Fzy Tzx count Fzx
    filename = sprintf('%s%s_HS.mat',directoryPath,myname);
    save(filename)
end

disp('done')

return
%% Ryx
% NOT READY
% fname = sprintf('%s_Ryx.bin',myname);
% fid = fopen(fname);
% [Data count] = fread(fid, H_mci.Ny*H_mci.Nx, 'float');
% fclose(fid);
% Ryx = reshape(Data,H_mci.Ny,H_mci.Nx);
% figure(5);clf
% imagesc(log10(Ryx))
