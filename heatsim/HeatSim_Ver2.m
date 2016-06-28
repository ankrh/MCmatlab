% Heat sim
% The purpose of this program is to simulate
% the heating of illuminated tissue.
% This program uses the light distribution and tissue structure from 
% maketissue_example2.m, makeTissueList.m, and lookmcxyz.m by importing
% HeatSimIn.mat

directoryPath='../Data/';

load([directoryPath 'Input_spectrum'])

%% Setup of simulation 
q = 0;
for fluence = F_min:dF:F_max
    q = q+1;
    num = num2str(q);
    load([directoryPath 'HeatSimIn_blood4_' num '.mat'])

%% Testing parameters


%% Simulation Parameters
%%% Parameters to set %%%
save_name           = ['HeatSimOut_blood4_' num '.mat'];
pulse_energy_area   = fluence; % [J/cm^2] pulse energy per area
pulse_duration      = pulse; % [s] pulse duration
duration_after      = 100e-3; % [s] simulation duration after pulse
dt                  = 5e-5; % timestep size, should be on the order of or less,
                            %than the smallest value of dx^2*HC/TC (characteristic timescale for heat diffusion in individual voxels
                            % also the pulse duration should be divisible by dt
Temp_initial        = 36; % initial temperature of tissue  
image_interval      = 1e-3; %[s] profiles will be saved with this interval

% swithches
model       = 4; % # of the tissue model loaded
save_on     = 1; % 1 to save output, 0 otherwise
Apcal       = 1; % 1 to simulate the time during illumination, 0 otherwise
postcal     = 1; % 1 to simulate the time after illumination, 0 otherwise
Bad_Temp    = 2000; % if the simulation reaches a temperature above this, something is wrong and the simulation terminates
quick_cal   = 1; % quick_cal calculates only on a smaller matrix, go to Quick Calculation Setup and tailor it to a 
                    % specific model

%%% Parameters that is set automatically %%%
cd ..
tissueProps = makeTissueList(nm); %Properties of the tissue types
cd heatsim
Temp        = Temp_initial*ones(400,400,400); % initial temperature in Celsius
[~,~,z_mesh]= meshgrid(1:400,1:400,0:dz:(Nz-1)*dz); % Temp gradient
Temp        = 30+6/((Nz-1)*dz)*z_mesh;              % for testing

% Matrices contaning the thermal properties of the tissue
HC = zeros(size(T));
TC = zeros(size(T));
for tissueNumber=1:length(tissueProps)
    HC(T==tissueNumber) = tissueProps(tissueNumber).HC;
    TC(T==tissueNumber) = tissueProps(tissueNumber).TC;
end

% miscellaneous constants
area            = (Nx*dx)^2; % surface area of the model
pulse_energy    = pulse_energy_area*area; % [J] energy delivered to the volume
Wdel            = pulse_energy/pulse_duration; % Watt delivered
Nt_light        = round(pulse_duration/dt); % Number of timesteps with illumination
Nt_no_light     = round(duration_after/dt); % Numer of timesteps with no illumination 
image_count     = 1; % a counter used for the 

%Due to the way photons leaving the matrix is handled in mcxyz.c, the Ap
%values on the borders are much too large, therefore they will be set to
%zero
Ap(1:Nx,1:Ny,1)     = 0;
Ap(1:Nx,1:Ny,Nz)    = 0;
Ap(1:Nx,1,1:Nz)     = 0;
Ap(1:Nx,Ny,1:Nz)    = 0;
Ap(1,1:Ny,1:Nz)     = 0;
Ap(Nx,1:Ny,1:Nz)    = 0;


%% Quick Calculation Setup
if quick_cal==1;
    if model ==3
        Nmax    = Nx;
        Nx      = 100;
        Ny      = 100;
        Nz      = 300;
        Nx_low  = Nmax/2+1-Nx/2;
        Nz_low  = 1;
        Nx_high = Nmax/2+Nx/2;
        Nz_high = Nz;
        z       = z(Nz_low:Nz_high);
        y       = y(Nx_low:Nx_high);
        x       = x(Nx_low:Nx_high);
        Temp    = Temp(Nx_low:Nx_high,Nx_low:Nx_high,Nz_low:Nz_high);
        Ap      = Ap(Nx_low:Nx_high,Nx_low:Nx_high,Nz_low:Nz_high);
        T       = T(Nx_low:Nx_high,Nx_low:Nx_high,Nz_low:Nz_high);
        HC      = HC(Nx_low:Nx_high,Nx_low:Nx_high,Nz_low:Nz_high);
        TC      = TC(Nx_low:Nx_high,Nx_low:Nx_high,Nz_low:Nz_high);
    elseif model==4
        Nmax    = Nx;
        Nx      = 100;
        Ny      = 100;
        Nz      = 100;
        Nx_low  = Nmax/2+1-Nx/2;
        Nz_low  = 1;
        Nx_high = Nmax/2+Nx/2;
        Nz_high = Nz;
        z       = z(Nz_low:Nz_high);
        y       = y(Nx_low:Nx_high);
        x       = x(Nx_low:Nx_high);
        Temp    = Temp(Nx_low:Nx_high,Nx_low:Nx_high,Nz_low:Nz_high);
        Ap      = Ap(Nx_low:Nx_high,Nx_low:Nx_high,Nz_low:Nz_high);
        T       = T(Nx_low:Nx_high,Nx_low:Nx_high,Nz_low:Nz_high);
        HC      = HC(Nx_low:Nx_high,Nx_low:Nx_high,Nz_low:Nz_high);
        TC      = TC(Nx_low:Nx_high,Nx_low:Nx_high,Nz_low:Nz_high);
    end
end

%% Heat Transfer Simulation with Illumination
if Apcal==1;
    tic
    for nt = 1:Nt_light
        % Calculates heat propagation
        dQx = (dt/dx)*dy*dz*((TC([2:Nx Nx],[1:Ny],[1:Nz])+TC([1:Nx],[1:Ny],[1:Nz]))./2.*(Temp([2:Nx Nx],[1:Ny],[1:Nz])-Temp([1:Nx],[1:Ny],[1:Nz]))+...
            -(TC([1:Nx],[1:Ny],[1:Nz])+TC([1 1:Nx-1],[1:Ny],[1:Nz]))./2.*(Temp([1:Nx],[1:Ny],[1:Nz])-Temp([1 1:Nx-1],[1:Ny],[1:Nz])));
        dQy =(dt/dy)*dx*dz*((TC([1:Nx],[2:Ny Ny],[1:Nz])+TC([1:Nx],[1:Ny],[1:Nz]))./2.*(Temp([1:Nx],[2:Ny Ny],[1:Nz])-Temp([1:Nx],[1:Ny],[1:Nz]))+...
            -(TC([1:Nx],[1:Ny],[1:Nz])+TC([1:Nx],[1 1:Ny-1],[1:Nz]))./2.*(Temp([1:Nx],[1:Ny],[1:Nz])-Temp([1:Nx],[1 1:Ny-1],[1:Nz])));
        dQz =(dt/dz)*dx*dy*((TC([1:Nx],[1:Ny],[2:Nz Nz])+TC([1:Nx],[1:Ny],[1:Nz]))./2.*(Temp([1:Nx],[1:Ny],[2:Nz Nz])-Temp([1:Nx],[1:Ny],[1:Nz]))+...
            -(TC([1:Nx],[1:Ny],[1:Nz])+TC([1:Nx],[1:Ny],[1 1:Nz-1]))./2.*(Temp([1:Nx],[1:Ny],[1:Nz])-Temp([1:Nx],[1:Ny],[1 1:Nz-1])));
        % Sum of heat propagation and heat generated from absorbed light
        dQ = dQx+dQy+dQz+Ap*dx*dy*dz*Wdel*dt; 
        % Calculating temperature at the next timestep
        Temp([1:Nx],[1:Ny],[2:Nz-1]) = Temp([1:Nx],[1:Ny],[2:Nz-1])+dQ([1:Nx],[1:Ny],[2:Nz-1])./HC([1:Nx],[1:Ny],[2:Nz-1])./(dx*dy*dz);
        % Break if the temperature increases too much
        if max(max(max(Temp)))>Bad_Temp
            break
        end
        % Saves temperature profiles
        if mod(dt*nt,image_interval)==0
            Tempzy(:,:,image_count)  = reshape(Temp(:,Nx/2,:),Ny,Nz)';
            image_count=image_count+1;
        end
        disp(['time: ' num2str(nt*dt) ', total time: ' num2str(Nt_light*dt)])
    end
    toc
end
Temp_post_light=Temp;
image_count_light_end=image_count-1;
%% Heat Transfer Simulation with No Illumination
if postcal==1
    tic
    for nt = 1:Nt_no_light
        dQx = (dt/dx)*dy*dz*((TC([2:Nx Nx],[1:Ny],[1:Nz])+TC([1:Nx],[1:Ny],[1:Nz]))./2.*(Temp([2:Nx Nx],[1:Ny],[1:Nz])-Temp([1:Nx],[1:Ny],[1:Nz]))+...
            -(TC([1:Nx],[1:Ny],[1:Nz])+TC([1 1:Nx-1],[1:Ny],[1:Nz]))./2.*(Temp([1:Nx],[1:Ny],[1:Nz])-Temp([1 1:Nx-1],[1:Ny],[1:Nz])));
        dQy =(dt/dy)*dx*dz*((TC([1:Nx],[2:Ny Ny],[1:Nz])+TC([1:Nx],[1:Ny],[1:Nz]))./2.*(Temp([1:Nx],[2:Ny Ny],[1:Nz])-Temp([1:Nx],[1:Ny],[1:Nz]))+...
            -(TC([1:Nx],[1:Ny],[1:Nz])+TC([1:Nx],[1 1:Ny-1],[1:Nz]))./2.*(Temp([1:Nx],[1:Ny],[1:Nz])-Temp([1:Nx],[1 1:Ny-1],[1:Nz])));
        dQz =(dt/dz)*dx*dy*((TC([1:Nx],[1:Ny],[2:Nz Nz])+TC([1:Nx],[1:Ny],[1:Nz]))./2.*(Temp([1:Nx],[1:Ny],[2:Nz Nz])-Temp([1:Nx],[1:Ny],[1:Nz]))+...
            -(TC([1:Nx],[1:Ny],[1:Nz])+TC([1:Nx],[1:Ny],[1 1:Nz-1]))./2.*(Temp([1:Nx],[1:Ny],[1:Nz])-Temp([1:Nx],[1:Ny],[1 1:Nz-1])));
        dQ=dQx+dQy+dQz;
        Temp([1:Nx],[1:Ny],[2:Nz-1])=Temp([1:Nx],[1:Ny],[2:Nz-1])+dQ([1:Nx],[1:Ny],[2:Nz-1])./HC([1:Nx],[1:Ny],[2:Nz-1])./(dx*dy*dz);
        % Saves temperature profiles
        if mod(dt*nt,image_interval)==0
            Tempzy(:,:,image_count)  = reshape(Temp(:,Nx/2,:),Ny,Nz)';
            image_count=image_count+1;
        end
        disp(['time: ' num2str(nt*dt) ', total time: ' num2str(Nt_no_light*dt)])
    end
    toc
end
%% Look at data
if model==4
   max(max(Temp_post_light([1:48 52:100],[1:48 52:100],13))) % highest temperature in epidermis 
end

maxTemp_post_light = max(max(max(Temp_post_light)))
%[maxTemp_papilla,I]=max(Tempzy(63,50,:)) %maximum temperature in the center of the papilla

Temp_post_lightzy = reshape(Temp_post_light(:,Nx/2,:),Ny,Nz)';
Tzx = reshape(T(Ny/2,:,:),Nx,Nz)';

% figure(1);clf
% imagesc(y,z,Tzx,[1 Nt])
% hold on
% cmap = makecmap(Nt);
% colormap(cmap)
% colorbar
% set(gca,'fontsize',18)
% set(colorbar,'fontsize',1)
% xlabel('x [cm]')
% ylabel('z [cm]')
% title('Tissue types')
% for i=1:Nt
%     yy = min(z) + (Nt-i)/(Nt-1)*max(z);
%     text(max(x)+1/6*(max(x)-min(x)),yy, sprintf('%d %s',i,tissue(i).s),'fontsize',12)
% end
% axis equal image
% 
% % draw launch
% N = 10; % # of beam rays drawn
% mcflag=0;
% switch mcflag
%     case 0 % uniform
%         for i=0:N
%             for j=-2:2
%             plot( [xs+radius*i/N xfocus + waist*j/2],[zs zfocus],'r-')
%             plot(-[xs+radius*i/N xfocus + waist*j/2],[zs zfocus],'r-')
%             end
%         end
% 
%     case 1 % Gaussian
%         for i=0:N
%             for j=-2:2
%             plot( [xs+radius*i/N xfocus + waist*j/2],[zs zfocus],'r-')
%             plot(-[xs+radius*i/N xfocus + waist*j/2],[zs zfocus],'r-')
%             end
%         end
% 
%     case 2 % iso-point
%         for i=1:20
%             th = (i-1)/19*2*pi;
%             xx = Nx/2*cos(th) + xs;
%             zz = Nx/2*sin(th) + zs;
%             plot([xs xx],[zs zz],'r-')
%         end
% end

% figure(5);clf
% imagesc(y,z,Temp_post_lightzy)
% hold on
% text(max(x)*0.9,min(z)-0.04*max(z),'T [^{\circ}C]','fontsize',18)
% colorbar
% set(gca,'fontsize',18)
% xlabel('y [cm]')
% ylabel('z [cm]')
% title('Temperature  [^{\circ}C] ')
% %colormap(makec2f)
% axis equal image
% 
% figure(6);clf
% imagesc(y,z,Tempzy(:,:,30))
% hold on
% text(max(x)*0.9,min(z)-0.04*max(z),'T [^{\circ}C]','fontsize',18)
% colorbar
% set(gca,'fontsize',18)
% xlabel('y [cm]')
% ylabel('z [cm]')
% title('Temperature  [^{\circ}C] ')
% %colormap(makec2f)
% axis equal image
% 
% figure(7);clf
% imagesc(y,z,Temp_post_lightzy(:,:)-Tempzy(:,:,1))
% hold on
% text(max(x)*0.9,min(z)-0.04*max(z),'T [^{\circ}C]','fontsize',18)
% colorbar
% set(gca,'fontsize',18)
% xlabel('y [cm]')
% ylabel('z [cm]')
% title('Temperature  [^{\circ}C] ')
% %colormap(makec2f)
% axis equal image
% %print -djpeg -r300 'Fig_Azy.jpg'
% 
% drawnow
if save_on==1
    save([directoryPath save_name],'Temp_post_light','Temp','Tempzy','pulse_energy_area','pulse_duration','duration_after','Nx','Ny','Nz','Nt','dx','dy','dz','T',...
        'x','y','z','tissue')
end
end

N=2500;
s=zeros(N,1);
for a=1:N
s(a)=tan(a); %*sin(-a/10);
end
Fs=2200; %increase value to speed up the sound, decrease to slow it down
soundsc(s,Fs) 
