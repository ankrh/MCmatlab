function lookmcxyz

% lookmcxyz.m
%   Looks at myname_F.bin, created by gomcxyz 
%   where myname is the name of the run: myname_T.bin, myname_H.mci
%
% Displays
%   Tzx = end-view of tissue structure
%   Fzx = end-view of vessel @ source, ys = -0.15 cm
%   Fzy = side-view along length of vessel
%
% Saves
%   Fzy_data4.mat = Fzy y z zzs Fdet
%       Fzy(400,400,8) = 8 z,y images
%       Fdet(8,1) = signal [1/cm^2] @ detector fiber
%

%% USER CHOICES <---------- you must specify -----
directoryPath = './Data/';
myname = 'blood4_broad_532';
nm     = 532;
saveon_HeatSim = 1;

%% Load header file
H_mci = reportHmci(directoryPath,myname);

format compact

%% Load Fluence rate F(y,x,z) 
filename = sprintf('%s%s_F.bin',directoryPath,myname);
disp(['loading ' filename])
tic
    fid = fopen(filename, 'rb');
    data = fread(fid, H_mci.Ny*H_mci.Nx*H_mci.Nz, 'float');
    fclose(fid);
toc
F = reshape(data,H_mci.Ny,H_mci.Nx,H_mci.Nz); % F(y,x,z)

% Load tissue structure in voxels, T(y,x,z) 
filename = sprintf('%s%s_T.bin',directoryPath,myname);
disp(['loading ' filename])
tic
    fid = fopen(filename, 'rb');
    data = fread(fid, H_mci.Ny*H_mci.Nx*H_mci.Nz, 'uint8');
    fclose(fid);
toc
T = reshape(data,H_mci.Ny,H_mci.Nx,H_mci.Nz); % T(y,x,z)

clear data

%%
x = ((1:H_mci.Nx)-H_mci.Nx/2-1/2)*H_mci.dx;
y = ((1:H_mci.Ny)-H_mci.Ny/2-1/2)*H_mci.dx;
z = ((1:H_mci.Nz)-1/2)*H_mci.dz;
ux = 2:H_mci.Nx-1;
uy = 2:H_mci.Ny-1;
uz = 2:H_mci.Nz-1;
zmin = min(z);
zmax = max(z);
zdiff = zmax-zmin;
xmin = min(x);
xmax = max(x);
xdiff = xmax-xmin;

%% Look at structure, Tzx
Tzx = reshape(T(H_mci.Ny/2,:,:),H_mci.Nx,H_mci.Nz)';
tissueList = makeTissueList(nm);
H_mci.Nt = length(tissueList);

figure(1); clf
imagesc(x(ux),z(uz),Tzx(uz,ux),[1 H_mci.Nt])
hold on
cmap = makecmap(H_mci.Nt);
colormap(cmap)
colorbar
set(gca,'fontsize',18)
set(colorbar,'fontsize',1)
xlabel('x [cm]')
ylabel('z [cm]')
title('Tissue types')
for i=1:H_mci.Nt
    yy = zmin + (H_mci.Nt-i)/(H_mci.Nt-1)*zdiff;
    text(xmin + xdiff*1.13,yy, sprintf('%d %s',i,tissueList(i).name),'fontsize',12)
end
axis equal image

%% draw launch
N = 10; % # of beam rays drawn
switch H_mci.mcflag
    case 0 % uniform
        for i=0:N
            for j=-2:2
            plot( [H_mci.xs+H_mci.radius*i/N H_mci.xfocus + H_mci.waist*j/2],[H_mci.zs H_mci.zfocus],'r-')
            plot(-[H_mci.xs+H_mci.radius*i/N H_mci.xfocus + H_mci.waist*j/2],[H_mci.zs H_mci.zfocus],'r-')
            end
        end

    case 1 % uniform over entire surface at height H_mci.zs
        for i=0:N
            plot( [xmin + (xmax-xmin)*i/N xmin + (xmax-xmin)*i/N],[H_mci.zs H_mci.zfocus],'r-')
        end

    case 2 % iso-point
        for i=1:20
            th = (i-1)/19*2*pi;
            xx = H_mci.Nx/2*cos(th) + H_mci.xs;
            zz = H_mci.Nx/2*sin(th) + H_mci.zs;
            plot([H_mci.xs xx],[H_mci.zs zz],'r-')
        end
end

name = sprintf('%s%s_tissue.jpg',directoryPath,myname);
print('-djpeg','-r300',name)


%% Look at Fluence Fzx @ launch point
Fzx = reshape(F(H_mci.Ny/2,:,:),H_mci.Nx,H_mci.Nz)'; % in z,x plane through source

figure(2);clf
imagesc(x,z,log10(Fzx),[-3 3])
hold on
text(max(x)*0.9,min(z)-0.04*max(z),'log_{10}( \phi )','fontsize',18)
colorbar
set(gca,'fontsize',18)
xlabel('x [cm]')
ylabel('z [cm]')
title('Fluence \phi [W/cm^2/W.delivered] ')
colormap(makec2f)
axis equal image
%axis([min(x) max(x) min(z) max(z)])

name = sprintf('%s%s_Fzx.jpg',directoryPath,myname);
print('-djpeg','-r300',name)

%% look Fzy
Fzy = reshape(F(:,H_mci.Nx/2,:),H_mci.Ny,H_mci.Nz)';

iy = round((H_mci.dy*H_mci.Ny/2 + 0.15)/H_mci.dy);
iz = round(H_mci.zs/H_mci.dz);
zzs  = H_mci.zs;
%Fdet = mean(reshape(Fzy(iz+[-1:1],iy+[0 1]),6,1));

figure(3);clf
imagesc(y,z,log10(Fzy),[-1 1]*3)
hold on
text(max(x)*0.9,min(z)-0.04*max(z),'log_{10}( \phi )','fontsize',18)
colorbar
set(gca,'fontsize',18)
xlabel('y [cm]')
ylabel('z [cm]')
title('Fluence \phi [W/cm^2/W.delivered] ')
colormap(makec2f)
axis equal image

name = sprintf('%s%s_Fzy.jpg',directoryPath,myname);
print('-djpeg','-r300',name)

drawnow
%% calculate Power Absorbtion

% F, matrix with fluency / input power
% T, Matrix containing tissue types
% tissueList, list of tissue properties tissueProps( tissue type, #) # = 1 is mua, # = 2 is mus, # = 3 is g

Ap = zeros(size(T));
for tissueNumber=1:length(tissueList)
   Ap(T==tissueNumber) = tissueList(tissueNumber).mua;
end
Azy = reshape(Ap(:,H_mci.Nx/2,:),H_mci.Ny,H_mci.Nz)';

figure(4);clf
imagesc(y,z,log10(Azy),[-1 1]*3)
hold on
text(max(x)*0.9,min(z)-0.04*max(z),'log_{10}( \phi )','fontsize',18)
colorbar
set(gca,'fontsize',18)
xlabel('y [cm]')
ylabel('z [cm]')
title('Power Absorbtion \phi [W/cm^3/W.delivered] ')
%colormap(makec2f)
axis equal image

name = sprintf('%s%s_Azy.jpg',directoryPath,myname);
print('-djpeg','-r300',name)

drawnow

% test H_mci.dx*H_mci.dy*H_mci.dz*sum(sum(sum(Ap)))<=1, to make sure that less than 1W is
% absorbed per 1W of incident power, this seems to be the case

%% Save HeatSim input

if saveon_HeatSim==1
    clearvars F Azy Fzy Tzx count Fzx
    filename = sprintf('%s%s_HS.mat',directoryPath,myname);
    save(filename)
end   

disp('done')

%% Ryx
% NOT READY
% fname = sprintf('%s_Ryx.bin',myname);
% fid = fopen(fname);
% [Data count] = fread(fid, H_mci.Ny*H_mci.Nx, 'float');
% fclose(fid);
% Ryx = reshape(Data,H_mci.Ny,H_mci.Nx);
% figure(5);clf
% imagesc(log10(Ryx))
