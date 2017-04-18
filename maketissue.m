function makeTissue

% maketissue.m
%   Creates a cube of optical property pointers,T(x,y,z), saved in
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
                        
directoryPath = 'Data/';
myname      = ['dentin_sim_' num2str(nm)];% name for files: myname_T.bin, myname_H.mci  
wavelength  = 850;      % [nm] set the range of wavelengths of the monte carlo simulation
time_min    = 0.5;      	% time duration of the simulation [min]
Nx = 100;               % # of bins in the x direction
Ny = Nx;                % # of bins in the y direction
Nz = Nx;                % # of bins in the z direction
dx = 1.0/Nx;            % size of x bins [cm]
dy = dx;                % size of y bins [cm]
dz = dx;                % size of z bins [cm]

% Set Monte Carlo launch flags
beamtypeflag = 5;     	% beam type: 0 = top-hat focus, top-hat far field beam,
                        % 1 = Gaussian focus, Gaussian far field beam,
                        % 2 = isotropically emitting point, 3 = infinite
                        % plane wave, 4 = pencil beam, 5 = top-hat focus,
                        % Gaussian far field beam, 6 = Gaussian focus,
                        % top-hat far field beam
boundaryflag = 1;       % 0 = no boundaries, 1 = escape at boundaries
                        % 2 = escape at surface only. No x, y, bottom z
                        % boundaries

% Set position of focus, only used for beamtypeflag ~=3 (if beamtypeflag == 2 this is the source position)
xfocus      = 0;        % set x position of focus
yfocus      = 0;        % set y position of focus
zfocus      = 1;    	% set z position of focus

% Set direction of beam center axis, only used if beamtypeflag ~= 2:
ux0         = 0;        % trajectory projected onto x axis
uy0         = 0;        % trajectory projected onto y axis
uz0         = sqrt(1 - ux0^2 - uy0^2); % such that ux^2 + uy^2 + uz^2 = 1, positive direction is assumed

% Set focus properties and divergence angles, only used if beamtypeflag == 0 or 1
waist       = 0.025;     % focus waist 1/e^2 radius in cm
divergence  = pi/8;    % divergence 1/e^2 half-angle of beam in rad
% divergence  = wavelength*1e-9/(pi*waist*1e-2); % Diffraction limited divergence angle for Gaussian beam

%% Prepare Monte Carlo 
format compact

% Create tissue properties
tissueList = makeTissueList(wavelength);
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

%% CREATE TISSUE STRUCTURE T(x,y,z)
%   Create T(x,y,z) by specifying a tissue type (an integer)
%   for each voxel in T.
%
%   Note: one need not use every tissue type in the tissue list.
%   The tissue list is a library of possible tissue types.

T = uint8(3*ones(Nx,Ny,Nz)); % fill background with air


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
%Enamel @ xc, zc, radius, oriented along x axis
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


%Dentin @ xc, zc, radius, oriented along x axis
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
     
%Blood @ xc, zc, radius, oriented along x axis
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

T = uint8(6*ones(Nx,Ny,Nz)); % fill with the testabsorber

%% Write the files
if SAVEON

    v = reshape(T,Nx*Ny*Nz,1);

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
        fprintf(fid,'%d\n'   ,beamtypeflag);
        fprintf(fid,'%d\n'   ,boundaryflag);
        fprintf(fid,'%0.8f\n',xfocus); % position of focus in cm
        fprintf(fid,'%0.8f\n',yfocus);
        fprintf(fid,'%0.8f\n',zfocus);
        fprintf(fid,'%0.8f\n',ux0); % beam center axis direction
        fprintf(fid,'%0.8f\n',uy0);
        fprintf(fid,'%0.8f\n',uz0);
        fprintf(fid,'%0.8f\n',waist); % waist radius in cm
        fprintf(fid,'%0.8f\n',divergence); % divergence half-angle of incoming beam in rad
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

figure(1);clf;
plotVolumetric(x,y,z,T,tissueList)
title('Tissue type illustration');

disp('done')

