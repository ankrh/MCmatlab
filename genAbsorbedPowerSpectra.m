function genAbsorbedPowerSpectra

% genAbsorbedPowerSpectra.m
%   loads at several myname_F.bin, created by 3Dmc.exe 
%   where myname is the name of the run: myname_T.bin, myname_H.mci
%   and myname end in _xxx where xxx=wavelength in nanometers
%   
%   Reads 'Elipse_12,2Jcm2_Pulse15ms_IntTime10s_OD3_mean' or similar file
%   to get measured data for light source spectrum
%
% For each run:
%   alc#_H.mci --> header:
%       timin,H_mci.Nx,H_mci.Ny,H_mci.Nz,dy,H_mci.dx,H_mci.dz,xs,ys,zx,H_mci.Nt,muav(),musv(),gv()
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

directoryPath = 'C:\Users\Kira Schmidt\Desktop\mcxyz';

load([directoryPath 'Input_spectrum'])
nm_s = nm;
P_s = Power;

m = 0;
for k = F_min:dF:F_max

    m = m+1;
    num = num2str(m);
    savename = ['HeatSimIn_dentin_' num '.mat'];


    %% Remember to define the size of Ap_total %%
    Ap_total=zeros(400,400,400);
    %%
    cc = 'rbgm'; % color
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

    for i_nm = length(Wavelengths):-1:1
        [~,I1] = min(abs(nm_s-(Wavelengths(i_nm)-dnm/2)));
        [~,I2] = min(abs(nm_s-(Wavelengths(i_nm)+dnm/2)));
        weight(i_nm) = trapz(nm_s(I1:I2),P_s(I1:I2));
    end
    weight = weight./sum(weight); % Normalise

    %% USER CHOICES <---------- you must specify -----
    for i_nm = 1:length(Wavelengths)
        nm = Wavelengths(i_nm);
        myname = ['dentin_sim_' num2str(nm)];

        %% Load header file
        H_mci = reportHmci(directoryPath,myname);

        %% Load Fluence rate F(y,x,z) 
        filename = sprintf('%s%s_F.bin',directoryPath,myname);
        disp(['loading ' filename])
        tic
            fid = fopen(filename, 'rb');
            [Data count] = fread(fid, H_mci.Ny*H_mci.Nx*H_mci.Nz, 'float');
            fclose(fid);
        toc
        F = reshape(Data,H_mci.Ny,H_mci.Nx,H_mci.Nz); % F(y,x,z)

        % Load tissue structure in voxels, T(y,x,z) 
        filename = sprintf('%s%s_T.bin',directoryPath,myname);
        disp(['loading ' filename])
        tic
            fid = fopen(filename, 'rb');
            [Data count] = fread(fid, H_mci.Ny*H_mci.Nx*H_mci.Nz, 'uint8=>uint8');
            fclose(fid);
        toc
        T = reshape(Data,H_mci.Ny,H_mci.Nx,H_mci.Nz); % T(y,x,z)

        clear Data

        %%
        x = ([1:H_mci.Nx]-H_mci.Nx/2-1/2)*H_mci.dx;
        y = ([1:H_mci.Ny]-H_mci.Ny/2-1/2)*H_mci.dx;
        z = ([1:H_mci.Nz]-1/2)*H_mci.dz;
        ux = [2:H_mci.Nx-1];
        uy = [2:H_mci.Ny-1];
        uz = [2:H_mci.Nz-1];
        zmin = min(z);
        zmax = max(z);
        zdiff = zmax-zmin;
        xmin = min(x);
        xmax = max(x);
        xdiff = xmax-xmin;
        tissueList = makeTissueList(nm);


        %% calculate Power Absorbtion

        % F, matrix with fluency / input power
        % T, Matrix containing tissue types
        % tissueList, list of tissue properties

        Ap = zeros(size(T));
        for tissueNumber=1:length(tissueList)
            Ap(T==tissueNumber) = tissueList(tissueNumber).mua;
        end
        Ap = Ap.*F;

        Ap_total = Ap_total+weight(i_nm)*Ap;
    end

    %% Save HeatSim input
    Ap=Ap_total;
    clearvars F Azy Fzy Tzx count Fzx Ap_total
    save([directoryPath savename])

end
