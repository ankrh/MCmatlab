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
format compact

%% USER CHOICES <---------- you must specify -----
directoryPath = 'exec/';
myname = 'dentin_sim_850';
wavelength = 850;
saveon_HeatSim = 1;

%% Load header file
H_mci = reportHmci(directoryPath,myname);
dx = H_mci.dx;
dy = H_mci.dy;
dz = H_mci.dz;
nx = H_mci.nx;
ny = H_mci.ny;
nz = H_mci.nz;

%% Load Fluence rate F(x,y,z)
filename = sprintf('%s%s_F.bin',directoryPath,myname);
disp(['Loading ' filename])
fid = fopen(filename,'rb');
data = fread(fid,nx*ny*nz,'float');
fclose(fid);
F = reshape(data,nx,ny,nz); % F(x,y,z)

%% Load tissue structure in voxels, T(x,y,z)
filename = sprintf('%s%s_T.bin',directoryPath,myname);
disp(['Loading ' filename])
fid = fopen(filename,'rb');
data = fread(fid,nx*ny*nz,'uint8=>uint8');
fclose(fid);
T = reshape(data,nx,ny,nz); % T(x,y,z)
clear data

%% Voxel center positions
x  = ((0:nx-1)-(nx-1)/2)*dx;
y  = ((0:ny-1)-(ny-1)/2)*dy;
z  = ((0:nz-1)+1/2)*dz;
tissueList = makeTissueList(wavelength);

%% Make volumetric tissue plot
figure(1);
plotVolumetric(x,y,z,T,tissueList);
title('Tissue type illustration');

%% Make volumetric fluence rate plot
figure(2);
plotVolumetric(x,y,z,F);
title('Fluence rate (Intensity) [W/cm^2/W.delivered] ')

%% Calculate normalized volumetric power (Power absorbed per cubic centimeter per watt delivered)
NVP = zeros(size(T));
for tissueNumber=1:length(tissueList)
    NVP(T==tissueNumber) = tissueList(tissueNumber).mua;
end
NVP = NVP.*F;

%% Make volumetric power absorption plot
figure(3);
plotVolumetric(x,y,z,NVP);
title('Normalized absorbed power per unit volume [W/cm^3/W.delivered] ')

%% Save HeatSim input
if saveon_HeatSim==1
    filename = sprintf('%s%s_HS.mat',directoryPath,myname);
    disp(['Saving ' filename])
    save(filename,'dx','dy','dz','T','tissueList','NVP')
end
return





