% Gen_Ap_mcxyz.m
%   loads at several myname_F.bin, created by 3Dmc.exe 
%   where myname is the name of the run: myname_T.bin, myname_H.mci
%   and myname end in _xxx where xxx=wavelength in nanometers
%   
%   Reads 'Elipse_12,2Jcm2_Pulse15ms_IntTime10s_OD3_mean' or similar file
%   to get measured data for light source spectrum
%
% For each run:
%   alc#_H.mci --> header:
%       timin,Nx,Ny,Nz,dy,dx,dz,xs,ys,zx,Nt,muav(),musv(),gv()
%   alc#_F.bin --> F(y,x,z) = relative fluence rate [1/cm^2]
%   alc#_T.bin --> T(y,x,z) = tissue types
%
% Displays
%   Tzx = end-view of tissue structure
%   Fzx = end-view of vessel @ source, ys = -0.15 cm
%   Fzy = side-view along length of vessel
%
% Saves
%   'HeatSimIn_model4.mat'
%

format compact

directoryPath = './Data/';

load([directoryPath 'Input_spectrum'])
nm_s = nm;
P_s = Power;

m = 0;
for k = F_min:dF:F_max
    
    m = m+1;
    num = num2str(m);
savename = ['HeatSimIn_blood4_' num '.mat'];
% savename = ['HeatSimIn_blood4_broad.mat'];


%% Remember to define the size of Ap_total %%
Ap_total=zeros(400,400,400);
%%
cc = 'rbgm'; % color
dnm = Dnm;
Wavelengths = nm_min:dnm:nm_max;
%% Load Spectra
% ext = '.mat'; % Only find files with this extention
% name = ['test2'];
% %name = ['spec_abs' numb 'J']; % Find all files which names include this string
% AllFiles = dirrec(char(fileparts(pwd)),ext); % dirrec finds all files in the folder and it's sub folder with the given extention
% n = 0;
% for i=1:length(AllFiles)
%     if isempty(strfind(char(AllFiles(i)),name))==0
%         if isempty(strfind(char(ext),'.mat'))==0 % Matlab data files
%             n=n+1;
%             files(n)=AllFiles(i);
%             load(char(AllFiles(i)));
%             nm_s(:,n) = nm(:);
%             P_s(:,n) = Power;
% %         elseif isempty(strfind(char(ext),'.ISD'))==0 % Raw measurement from the sphere lab
% %             n=n+1;
% %             files(n)=AllFiles(i);
% %             Import = importdata(char(files(n)));
% %             data = Import.data;
% %             text(:,n) = Import.textdata;
% %             nm_s(:,n) = data(:,1);
% %             P_s(:,n) = data(:,2);
%         end
%     end
% end
% clear nm Power Import data
% if n>1; display('The name specified for the file to load the spectra from was not unique'); break; end
% if max(Wavelengths+dnm/2) > max(nm_s) || min(Wavelengths-dnm/2) < min(nm_s)
%     display('The spectra doesn''t cover all the wavelengths specified');break; 
% end

for i = 1:length(Wavelengths) 
    [~,I1] = min(abs(nm_s-(Wavelengths(i)-dnm/2)));
    [~,I2] = min(abs(nm_s-(Wavelengths(i)+dnm/2)));
    weight(i) = trapz(nm_s(I1:I2),P_s(I1:I2));
end
weight = weight./sum(weight); % Normalise

%%%% USER CHOICES <---------- you must specify -----
index=1;
%weight = ones(1,6)/6;
for i = 1:length(Wavelengths)
    nm = Wavelengths(i);
myname = ['blood4_broad_' num2str(nm)];
%%%%


disp(sprintf('------ mcxyz %s -------',myname))

% Load header file
filename = sprintf('%s%s_H.mci',directoryPath,myname);
disp(['loading ' filename])
fid = fopen(filename, 'r');
A = fscanf(fid,'%f',[1 Inf])';
fclose(fid);

%% parameters
time_min = A(1);
Nx = A(2);
Ny = A(3);
Nz = A(4);
dx = A(5);
dy = A(6);
dz = A(7);
mcflag = A(8);
launchflag = A(9);
xs = A(10);
ys = A(11);
zs = A(12);
xfocus = A(13);
yfocus = A(14);
zfocus = A(15);
ux0 = A(16);
uy0 = A(17);
uz0 = A(18);
radius = A(19);
waist = A(20);
Nt = A(21);
j = 21;
for i=1:Nt
    j=j+1;
    muav(i,1) = A(j);
    j=j+1;
    musv(i,1) = A(j);
    j=j+1;
    gv(i,1) = A(j);
end

reportHmci(directoryPath,myname)

%% Load Fluence rate F(y,x,z) 
filename = sprintf('%s%s_F.bin',directoryPath,myname);
disp(['loading ' filename])
tic
    fid = fopen(filename, 'rb');
    [Data count] = fread(fid, Ny*Nx*Nz, 'float');
    fclose(fid);
toc
F = reshape(Data,Ny,Nx,Nz); % F(y,x,z)


% Load tissue structure in voxels, T(y,x,z) 
filename = sprintf('%s%s_T.bin',directoryPath,myname);
disp(['loading ' filename])
tic
    fid = fopen(filename, 'rb');
    [Data count] = fread(fid, Ny*Nx*Nz, 'uint8');
    fclose(fid);
toc
T = reshape(Data,Ny,Nx,Nz); % T(y,x,z)

clear Data

%%
x = ([1:Nx]-Nx/2-1/2)*dx;
y = ([1:Ny]-Ny/2-1/2)*dx;
z = ([1:Nz]-1/2)*dz;
ux = [2:Nx-1];
uy = [2:Ny-1];
uz = [2:Nz-1];
zmin = min(z);
zmax = max(z);
zdiff = zmax-zmin;
xmin = min(x);
xmax = max(x);
xdiff = xmax-xmin;
cd ..
tissueProps = makeTissueList(nm);
cd heatsim


%% calculate Power Absorbtion

% F, matrix with fluency / input power
% T, Matrix containing tissue types
% tissueProps, list of tissue properties tissueProps( tissue type, #) # = 1 is mua, # = 2 is mus, # = 3 is g
% 
for tissueNumber=1:length(tissueProps)
    T(T==tissueNumber) = tissueProps(tissueNumber).mua;
end
Ap = T.*F;

Ap_total = Ap_total+weight(index)*Ap;
index=index+1;
end

Azy  = reshape(Ap_total(:,Nx/2,:),Ny,Nz)';

% figure(4);clf
% imagesc(y,z,log10(Azy),[-1 1]*3)
% hold on
% text(max(x)*0.9,min(z)-0.04*max(z),'log_{10}( \phi )','fontsize',18)
% colorbar
% set(gca,'fontsize',18)
% xlabel('y [cm]')
% ylabel('z [cm]')
% title('Power Absorbtion \phi [W/cm^3/W.delivered] ')
% % colormap(makec2f)
% axis equal image
% 
% print -djpeg -r300 'Fig_Azy.jpg'
% 
% drawnow

% test dx*dy*dz*sum(sum(sum(Ap)))<=1, to make sure that less than 1W is
% absorbed per 1W of incident power, this seems to be the case

%% Save HeatSim input
Ap=Ap_total;
clearvars F Azy Fzy Tzx count Fzx Ap_total
save([directoryPath savename])

end

N=2500;
s=zeros(N,1);
for a=1:N
s(a)=tan(a); %*sin(-a/10);
end
Fs=2200; %increase value to speed up the sound, decrease to slow it down
soundsc(s,Fs) 
