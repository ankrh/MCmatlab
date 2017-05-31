function makeTissue

% maketissue.m
%   Creates a cube of optical property pointers,T(x,y,z), saved in
%       myname_T.bin = a tissue structure file
%   which specifies a complex tissue for use by mcxyz.c.
%
%   Also prepares a listing of the optical properties at chosen wavelength
%   for use by mcxyz, [mua, mus, g], for each tissue type specified
%   in myname_T.bin. This listing is saved in
%       myname_H.mci = the input file for use by mcxyz.
%
%   Will generate a figure illustrating the tissue with its various
%   tissue types.
%
%   Uses
%       makeTissueList.m
%
%   To use, 
%       1. Prepare makeTissueList.m so that it contains the tissue
%   types desired.
%       2. Specify the USER CHOICES below.
%       2. Run this program, maketissue.m.
%
%   Note: mcxyz uses optical properties in cm^-1.
%       

%% USER CHOICES <-------- You must set these parameters ------
SAVEON      = 1;        % 1 = save myname_T.bin, myname_H.mci 
                        % 0 = don't save. Just check the program.
                        
directoryPath = 'exec/';
wavelength  = 850;      % [nm] set the range of wavelengths of the monte carlo simulation
myname      = ['hair_sim_' num2str(wavelength)];% name for files: myname_T.bin, myname_H.mci  
simulationTimeRequested    = 0.5;      	% time duration of the simulation [min]
nx = 100;               % # of bins in the x direction
ny = nx;                % # of bins in the y direction
nz = nx;                % # of bins in the z direction
dx = 0.05/nx;            % size of x bins [cm]
dy = 0.05/ny;                % size of y bins [cm]
dz = 0.15/nz;                % size of z bins [cm]

% Set Monte Carlo launch flags
beamtypeFlag = 3;     	% beam type: 0 = top-hat focus, top-hat far field beam,
                        % 1 = Gaussian focus, Gaussian far field beam,
                        % 2 = isotropically emitting point, 3 = infinite
                        % plane wave, 4 = pencil beam, 5 = top-hat focus,
                        % Gaussian far field beam, 6 = Gaussian focus,
                        % top-hat far field beam
boundaryFlag = 1;       % 0 = no boundaries, 1 = escape at boundaries
                        % 2 = escape at surface only. No x, y, bottom z
                        % boundaries

% Set position of focus, only used for beamtypeflag ~=3 (if beamtypeflag == 2 this is the source position)
xFocus      = 0;        % set x position of focus
yFocus      = 0;        % set y position of focus
zFocus      = 1;    	% set z position of focus

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

% Specify Monte Carlo parameters    
x  = ((0:nx-1)-(nx-1)/2)*dx;
y  = ((0:ny-1)-(ny-1)/2)*dy;
z  = ((0:nz-1)+1/2)*dz;

if isinf(zFocus), zFocus = 1e12; end

%% CREATE TISSUE STRUCTURE T(x,y,z)
%   Create T(x,y,z) by specifying a tissue type (an integer)
%   for each voxel in T.
%
%   Note: one need not use every tissue type in the tissue list.
%   The tissue list is a library of possible tissue types.

T = 9*ones(nx,ny,nz,'uint8'); % fill background with air

zsurf = 0.02;  % position of gel/skin surface[cm]
epd_thick = 0.01; %Thickness of the epidermis [cm]

hair_diameter = 0.0075; % varies from 17 - 180 micrometers, should increase with colouring and age
hair_radius = hair_diameter/2;      	% hair radius [cm]
hair_bulb_radius = 1.7*hair_radius; % [cm]
ratio_papilla=5/12;
papilla_radius = hair_bulb_radius*ratio_papilla;
hair_depth = 0.1; % varies from 0.06-0.3cm
xi_start = nx/2-round(hair_radius/dx);
xi_end = nx/2+round(hair_radius/dx);
yi_start = ny/2-round(hair_radius/dy);
yi_end = ny/2+round(hair_radius/dy);

xi_start2 = nx/2-round(hair_bulb_radius/dx);
xi_end2 = nx/2+round(hair_bulb_radius/dx);
yi_start2 = ny/2-round(hair_bulb_radius/dy);
yi_end2 = ny/2+round(hair_bulb_radius/dy);

for iz=1:nz % for every depth z(iz)
    
    % epidermis
    if iz>round((zsurf)/dz) && iz<=round((zsurf+epd_thick)/dz)
        T(:,:,iz) = 10;
    end
    % Gel
    if iz<=round(zsurf/dz)
        T(:,:,iz) = 4;
    end
    
    % Hair @ xc, yc, radius, oriented along z axis
    xc      = 0;            % [cm], center of hair
    yc      = 0;     	% [cm], center of hair
    zc      = zsurf+hair_depth; % center of hair bulb
    zc_papilla = zsurf+hair_depth+sqrt(2)*(1-ratio_papilla)*hair_bulb_radius; %[cm] z-coordinate center of the papilla
    if iz<=round((zsurf+hair_depth)/dz)
        for ix=xi_start:xi_end
            for iy=yi_start:yi_end
                xd = x(ix) - xc;	% vessel, x distance from vessel center
                yd = y(iy) - yc;   	% vessel, z distance from vessel center
                zd = z(iz)-zc;
                r  = sqrt(xd^2 + yd^2);	% radius from vessel center
                if (r<=hair_radius)     	% if r is within hair
                    T(iy,ix,iz) = 15; % hair
                end
            end % iy
            
        end %ix
    end
    
    %Hair Bulb
    for ix=xi_start2:xi_end2
        for iy=yi_start2:yi_end2
            xd = x(ix) - xc;	% vessel, x distance from vessel center
            yd = y(iy) - yc;   	% vessel, z distance from vessel center
            zd = z(iz)-zc;
            r2 = sqrt(xd^2 + yd^2 + 1/2*zd^2); % radius from bulb center
            if (r2<=hair_bulb_radius)     	% if r2 is within hair bulb
                T(iy,ix,iz) = 15; % hair
            end
        end % iy
        
    end %ix
    
    %Papilla
    for ix=xi_start2:xi_end2
        for iy=yi_start2:yi_end2
            xd = x(ix) - xc;	% vessel, x distance from vessel center
            yd = y(iy) - yc;   	% vessel, z distance from vessel center
            zd = z(iz)-zc_papilla;
            r3 = sqrt(xd^2 + yd^2 + 1/2*zd^2); % radius from papilla center
            if (r3<=papilla_radius)     	% if r2 is within hair bulb
                T(iy,ix,iz) = 9; % dermis, standin for papilla tissue
            end
        end % iy
        
    end %ix
    
end %iz

%% Discard the unused tissue types
% usedTissues = unique(T);
% newTissueNumber = 0;
% for usedTissueNumber = usedTissues'
%     newTissueNumber = newTissueNumber + 1;
%     T(T(:)==usedTissueNumber) = newTissueNumber;
%     tissueList(newTissueNumber) = tissueList(usedTissueNumber);
% end
% tissueList = tissueList(1:newTissueNumber);

%% Write the files
if SAVEON

    Nt = length(tissueList);
    for i=Nt:-1:1
        muav(i)  = tissueList(i).mua;
        musv(i)  = tissueList(i).mus;
        gv(i)    = tissueList(i).g;
    end

    v = reshape(T,nx*ny*nz,1);

    %% WRITE FILES
    % Write myname_H.mci file
    %   which contains the Monte Carlo simulation parameters
    %   and specifies the tissue optical properties for each tissue type.
    commandwindow
    fprintf('--------create %s --------\n',myname)
    filename = sprintf('%s%s_H.mci',directoryPath,myname);
    fid = fopen(filename,'w');
        % run parameters
        fprintf(fid,'%0.2f\n',simulationTimeRequested);
        fprintf(fid,'%d\n'   ,nx);
        fprintf(fid,'%d\n'   ,ny);
        fprintf(fid,'%d\n'   ,nz);
        fprintf(fid,'%0.8f\n',dx);
        fprintf(fid,'%0.8f\n',dy);
        fprintf(fid,'%0.8f\n',dz);
        % launch parameters
        fprintf(fid,'%d\n'   ,beamtypeFlag);
        fprintf(fid,'%d\n'   ,boundaryFlag);
        fprintf(fid,'%0.8f\n',xFocus); % position of focus in cm
        fprintf(fid,'%0.8f\n',yFocus);
        fprintf(fid,'%0.8f\n',zFocus);
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

