function makeTissue

% maketissue.m
%   Creates a cube of optical property pointers,T(y,x,z), saved in
%       myname_T.bin = a tissue structure file
%   which specifies a complex tissue for use by mcxyz.c.
%
%   Also prepares a listing of the optical properties at chosen wavelength
%   for use by mcxyz.c, [mua, mus, g], for each tissue type specified
%   in myname_T.bin. This listing is saved in
%       myname_H.mci = the input file for use by mcxyz.c.
%
%   Will generate a figure illustrating the tissue with its various
%   tissue types and the beam being launched.
%
%   Uses
%       makeTissueList.m
%
%   To use, 
%       1. Prepare makeTissueList.m so that it contains the tissue
%   types desired.
%       2. Specify the USER CHOICES.
%       2. Run this program, maketissue.m.
%
%   Note: mcxyz.c can use optical properties in cm^-1 or mm^-1 or m^-1,
%       if the bin size (binsize) is specified in cm or mm or m,
%       respectively.
%
%  Steven L. Jacques. updated Aug 21, 2014.
%       

%% USER CHOICES <-------- You must set these parameters ------
SAVEON      = 1;        % 1 = save myname_T.bin, myname_H.mci 
                        % 0 = don't save. Just check the program.
                        
nm          = 850;       % set the range of wavelengths of the monte carlo simulation
directoryPath = 'Data/';
myname      = ['dentin_sim_' num2str(nm)];% name for files: myname_T.bin, myname_H.mci  
time_min    = 5;      	% time duration of the simulation [min]
Nx = 250;               % # of bins in the x direction
Ny = Nx;                % # of bins in the y direction
Nz = Nx;                % # of bins in the z direction
dx = 1.0/Nx;            % size of x bins [cm]
dy = dx;                % size of y bins [cm]
dz = dx;                % size of z bins [cm]

% Set Monte Carlo launch flags
mcflag      = 1;     	% launch: 0 = top hat beam, 1 = infinite plane wave over entire surface, 2 = isotropic pt. 
launchflag  = 0;        % 0 = let mcxyz.c calculate launch trajectory
                        % 1 = manually set launch vector.
boundaryflag = 1;       % 0 = no boundaries, 1 = escape at boundaries
                        % 2 = escape at surface only. No x, y, bottom z
                        % boundaries

% Sets position of source
xs          = 0;      	% x of source [cm]
ys          = 0;        % y of source [cm]
zs          = 0;        % z of source [cm]

% Set position of focus, so mcxyz can calculate launch trajectory
xfocus      = 0;        % set x,position of focus
yfocus      = 0;        % set y,position of focus
zfocus      = inf;    	% set z,position of focus (=inf for collimated beam)

% only used if mcflag == 0 (top hat beam)
radius      = 0.05;      % 1/e radius of beam at tissue surface
waist       = 0.010;  	% 1/e radius of beam at focus

% only used if launchflag == 1 (manually set launch trajectory):
ux0         = 0.7;      % trajectory projected onto x axis
uy0         = 0.4;      % trajectory projected onto y axis
uz0         = sqrt(1 - ux0^2 - uy0^2); % such that ux^2 + uy^2 + uz^2 = 1

%% Prepare Monte Carlo 
format compact

% Create tissue properties
tissueList = makeTissueList(nm);
Nt = length(tissueList);
for i=Nt:-1:1
    muav(i)  = tissueList(i).mua;
    musv(i)  = tissueList(i).mus;
    gv(i)    = tissueList(i).g;
end

% Specify Monte Carlo parameters    
x  = ((0:Nx-1)-(Nx-1)/2)*dx;
y  = ((0:Ny-1)-(Ny-1)/2)*dy;
z  = ((0:Nz-1)+1/2)*dz;
xmin = min(x);
xmax = max(x);

if isinf(zfocus), zfocus = 1e12; end

%% CREATE TISSUE STRUCTURE T(y,x,z)
%   Create T(y,x,z) by specifying a tissue type (an integer)
%   for each voxel in T.
%
%   Note: one need not use every tissue type in the tissue list.
%   The tissue list is a library of possible tissue types.

T = uint8(3*ones(Ny,Nx,Nz)); % fill background with air

zsurf       = 0.01;  % position of [cm]
dentin_depth = 0.26;
enamel_depth = 0.001; % Enamel depth [cm]
dentin_thick = 0.25; %Thickness of the dentin [cm]
enamel_thick = 0.25; %thickness of enamel [cm]

enamel_diameter = 0.2; % varies from 17 - 180 micrometers, should decrease with age
enamel_radius = enamel_diameter/2;      	% hair radius [cm]

xi_start = Nx/2-round(enamel_radius/dx);
xi_end = Nx/2+round(enamel_radius/dx);
yi_start = Ny/2-round(enamel_radius/dy);
yi_end = Ny/2+round(enamel_radius/dy);

xi_start2 = Nx/2-round(enamel_radius/dx);
xi_end2 = Nx/2+round(enamel_radius/dx);
yi_start2 = Ny/2-round(enamel_radius/dy);
yi_end2 = Ny/2+round(enamel_radius/dy);

for iz=1:Nz % for every depth z(iz)
 
    
% enamel and dentin for only cross sectional model     
    %enamel
    %if iz>round((zsurf+enamel_depth)/dz) && iz<=round((zsurf+enamel_depth+enamel_thick)/dz)
    %   T(:,:,iz) = 2;
    %end
    % dentin
    %if iz>round((zsurf+dentin_depth)/dz) && iz<=round((zsurf+dentin_depth+dentin_thick)/dz)
    %    T(:,:,iz) = 1;
    %end
    

    
%     
%Enamel @ xc, zc, radius, oriented along y axis
     yc      = 0;            % [cm], center of blood vessel
     zc      = Nz/2*dz;     	% [cm], center of blood vessel
    enamelradius  = 0.300;      	% blood vessel radius [cm]
     for iy=1:Ny
             yd = y(iy) - yc;	% vessel, x distance from vessel center
             zd = z(iz) - zc;   	% vessel, z distance from vessel center                
             r  = sqrt(yd^2 + zd^2);	% r from vessel center
             if (r<=enamelradius)     	% if r is within vessel
                 T(:,iy,iz) = 2; % blood
             end
 
     end %ix


%Dentin @ xc, zc, radius, oriented along y axis
     yc      = 0;            % [cm], center of blood vessel
     zc      = Nz/2*dz;     	% [cm], center of blood vessel
     dentinradius  = 0.200;      	% blood vessel radius [cm]
     for iy=1:Ny
             yd = y(iy) - yc;	% vessel, x distance from vessel center
             zd = z(iz) - zc;   	% vessel, z distance from vessel center                
             r  = sqrt(yd^2 + zd^2);	% r from vessel center
             if (r<=dentinradius)     	% if r is within vessel
                 T(:,iy,iz) = 1; % blood
             end
 
     end %ix
     
%Blood @ xc, zc, radius, oriented along y axis
     yc      = 0;            % [cm], center of blood vessel
     zc      = Nz/2*dz;     	% [cm], center of blood vessel
     vesselradius  = 0.050;      	% blood vessel radius [cm]
     for iy=1:Ny
             yd = y(iy) - yc;	% vessel, x distance from vessel center
             zd = z(iz) - zc;   	% vessel, z distance from vessel center                
             r  = sqrt(yd^2 + zd^2);	% r from vessel center
             if (r<=vesselradius)     	% if r is within vessel
                 T(:,iy,iz) = 5; % blood
             end
 
     end %ix     

    
% Hair @ xc, yc, radius, oriented along z axis
%      xc      = 0.02;            % [cm], center of hair
%      yc      = 0.02;     	% [cm], center of hair
%      zc      = zsurf+enamel_depth; % center of hair bulb
%      %zc_enamel = zsurf+enamel_depth+sqrt(2)*(1-ratio_enamel)*enamel_radius; %[cm] z-coordinate center of the papilla
%      %if iz>round(zsurf/dz) && iz<=round((zsurf+enamel_depth)/dz)
%          for ix=xi_start:xi_end
%             for iy=yi_start:yi_end
%                  xd = x(ix) - xc;	% vessel, x distance from vessel center
%                  yd = y(iy) - yc;   	% vessel, z distance from vessel center
%                  zd = z(iz)-zc;
%                  r  = sqrt(xd^2 + yd^2);	% radius from vessel center
%                  if (r<=enamel_radius)     	% if r is within hair
%                      T(iy,ix,iz) = 2; % enamel 
%                  end
%              end % iy
%              
%          end %ix
%      end  
    
%blood vessel @ xc, zc, radius, oriented along y axis
     %xc      = 0;            % [cm], center of blood vessel
     %zc      = Nz/2*dz;     	% [cm], center of blood vessel
     %vesselradius  = 0.100;      	% blood vessel radius [cm]
     %for ix=1:Nx
        %     xd = x(ix) - xc;	% vessel, x distance from vessel center
         %    zd = z(iz) - zc;   	% vessel, z distance from vessel center                
         %    r  = sqrt(xd^2 + zd^2);	% r from vessel center
          %   if (r<=vesselradius)     	% if r is within vessel
           %      T(:,ix,iz) = 1; % blood
           %  end
 
     %end %ix
%    
% 
%
%     % Hair Bulb
%    for ix=xi_start2:xi_end2
%        for iy=yi_start2:yi_end2
%              xc = 0.02; 
%              yc = 0.02;
%              zc = 0.01;
%              xd = x(ix) - xc;	% vessel, x distance from vessel center
%              yd = y(iy) - yc;   	% vessel, z distance from vessel center
%              zd = z(iz)-zc;
%              r2 = sqrt(xd^2 + yd^2 + 1/2*zd^2); % radius from bulb center
%              if (r2<=enamel_radius)     	% if r2 is within hair bulb
%                  T(iy,ix,iz) = 4; % hair
%              end
%          end % iy
%          
%      end %ix
%
%     % Papilla
%     for ix=xi_start2:xi_end2
%         for iy=yi_start2:yi_end2
%             xd = x(ix) - xc;	% vessel, x distance from vessel center
%             yd = y(iy) - yc;   	% vessel, z distance from vessel center
%             zd = z(iz)-zc_papilla;
%             r3 = sqrt(xd^2 + yd^2 + 1/2*zd^2); % radius from papilla center
%             if (r3<=papilla_radius)     	% if r2 is within hair bulb
%                 T(iy,ix,iz) = 4; % dermis, standin for papilla tissue
%             end
%         end % iy
%         
%     end %ix
end % iz

%T= shiftdim(T,1); % shifts dimension to have the light coming from side

%% Write the files
if SAVEON

    v = reshape(T,Ny*Nx*Nz,1);

    %% WRITE FILES
    % Write myname_H.mci file
    %   which contains the Monte Carlo simulation parameters
    %   and specifies the tissue optical properties for each tissue type.
    commandwindow
    fprintf('--------create %s --------\n',myname)
    filename = sprintf('%s%s_H.mci',directoryPath,myname);
    fid = fopen(filename,'w');
        % run parameters
        fprintf(fid,'%0.2f\n',time_min);
        fprintf(fid,'%d\n'   ,Nx);
        fprintf(fid,'%d\n'   ,Ny);
        fprintf(fid,'%d\n'   ,Nz);
        fprintf(fid,'%0.8f\n',dx);
        fprintf(fid,'%0.8f\n',dy);
        fprintf(fid,'%0.8f\n',dz);
        % launch parameters
        fprintf(fid,'%d\n'   ,mcflag);
        fprintf(fid,'%d\n'   ,launchflag);
        fprintf(fid,'%d\n'   ,boundaryflag);
        fprintf(fid,'%0.8f\n',xs);
        fprintf(fid,'%0.8f\n',ys);
        fprintf(fid,'%0.8f\n',zs);
        fprintf(fid,'%0.8f\n',xfocus);
        fprintf(fid,'%0.8f\n',yfocus);
        fprintf(fid,'%0.8f\n',zfocus);
        fprintf(fid,'%0.8f\n',ux0); % if manually setting ux,uy,uz
        fprintf(fid,'%0.8f\n',uy0);
        fprintf(fid,'%0.8f\n',uz0);
        fprintf(fid,'%0.8f\n',radius);
        fprintf(fid,'%0.8f\n',waist);
        % tissue optical properties
        fprintf(fid,'%d\n',Nt);
        for i=1:Nt
            fprintf(fid,'%0.8f\n',muav(i));
            fprintf(fid,'%0.8f\n',musv(i));
            fprintf(fid,'%0.8f\n',gv(i));
        end
    fclose(fid);

    %% write myname_T.bin file
    filename = sprintf('%s%s_T.bin',directoryPath,myname);
    disp(['create ' filename])
    fid = fopen(filename,'wb');
    fwrite(fid,v,'uint8');
    fclose(fid);

end % SAVEON

%% Look at structure of Tzx at iy=Ny/2
Tzx  = squeeze(T(Ny/2,:,:))'; % Tyxz -> Txz -> Tzx

figure(1); clf
plotTissue(Tzx,tissueList,x,z)
hold on

%% draw launch
switch mcflag
    case 0 % uniform
        for i=0:5
            for j=-2:2
                plot( [xs+radius*i/5 xfocus + waist*j/2],[zs zfocus],'r-')
                plot(-[xs+radius*i/5 xfocus + waist*j/2],[zs zfocus],'r-')
            end
        end

    case 1 % uniform over entire surface at height zs
        for i=0:10
            plot( [xmin + (xmax-xmin)*i/10 xmin + (xmax-xmin)*i/10],[zs zfocus],'r-')
        end

    case 2 % iso-point
        for i=1:20
            th = (i-1)/19*2*pi;
            xx = Nx/2*cos(th) + xs;
            zz = Nx/2*sin(th) + zs;
            plot([xs xx],[zs zz],'r-')
        end
end

axis([min(x) max(x) min(z) max(z)])

figure(2);clf;
plotVolumetric(x,y,z,T,tissueList)
title('Tissue type illustration');

disp('done')

