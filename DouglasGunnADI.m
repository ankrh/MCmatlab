%% User specified parameters
calctype = 4; % 1: 1D, 2:2D, 3:3D, 4:mex 3D

heatsinkedboundary = false;

dt = 0.0003;

TC = [0.0037];

VHC = [3.7606];

nx = 100;
ny = 80;
nz = 60;

dx = 0.1/nx;
dy = 0.11/ny;
dz = 0.12/nz;

nt = 10;

%% Calculations
nM = length(TC);
switch calctype
    case 1
		%% Initialize arrays
		x = dx*(-(nx-1)/2:(nx-1)/2);
        N = nx;
        M   = ones(nx,1);
        M(12:13) = 2;
        ABS = zeros(nx,1); % Power deposited per unit length/area/volume per unit time
%         ABS(10) = 100;
		
		A1 = NaN(nM,nM,1);
		s = zeros(N,1);
		for i=1:nM
			for j=1:nM
				TC_eff = 2*TC(i)*TC(j)/(TC(i)+TC(j));
				if isnan(TC_eff); TC_eff = 0; end
				A1(i,j) = dt*TC_eff/(2*VHC(i)*dx^2); % A1 coupling factor for heat flow from medium j into medium i
			end
		end
		for i=1:N
			ix = i;
            heatsinkedvoxel = heatsinkedboundary && (ix == 1 || ix == nx);
			if ~heatsinkedvoxel
				s(i) = ABS(i)*dt/VHC(M(i));
			end
		end

		%% Initialize figure
        T = zeros(N,1);
		d = zeros(N,1);
		T(20) = 1;
		b = NaN(nx,1);
		figure(1);clf;
        h_line = plot(T);
        h_title = title('0');
		
		%% Do loop
        for it=1:nt
			%% Explicit part of the step:
			for ix = 1:nx % construct d, the right hand side of the system of linear equations
				if ~(heatsinkedboundary && (ix == 1 || ix == nx))
					delta = s(ix);
					if ix ~= 1;  delta = delta + A1(M(ix),M(ix-1))*(T(ix-1)-T(ix)); end
					if ix ~= nx; delta = delta + A1(M(ix),M(ix+1))*(T(ix+1)-T(ix)); end
				else
					delta = 0;
				end
				d(ix) = T(ix) + delta;
			end
			
			%% Implicit part of the step, Thomas algorithm, sweeps up from 2 to nx and then down from nx to 1
			% Forward sweep, index 2:
			if heatsinkedboundary
				b(1) = 1;
			else
				b(1) = 1 + A1(M(1),M(2));
			end
			w = -A1(M(2),M(1))/b(1);
			if heatsinkedboundary
				b(2) = 1 + A1(M(2),M(3)) + A1(M(2),M(1));
			else
				b(2) = 1 + A1(M(2),M(3)) + A1(M(2),M(1)) + w*A1(M(1),M(2));
			end
			d(2) = d(2) - w*d(1);
			% Forward sweep, indices 3 to nx-1:
			for ix = 3:nx-1
				w = -A1(M(ix),M(ix-1))/b(ix-1);
				b(ix) = 1 + A1(M(ix),M(ix-1)) + A1(M(ix),M(ix+1));
				b(ix) = b(ix) + w*A1(M(ix-1),M(ix));
				d(ix) = d(ix) - w*d(ix-1);
			end
			% Forward sweep, index nx:
			if heatsinkedboundary
				b(nx) = 1;
			else
				w = -A1(M(nx),M(nx-1))/b(nx-1);
				b(nx) = 1 + A1(M(nx),M(nx-1));
				b(nx) = b(nx) + w*A1(M(nx-1),M(nx));
				d(nx) = d(nx) - w*d(nx-1);
			end

			% Back sweep, index nx:
			T(nx) = d(nx)/b(nx);
			% Back sweep, indices nx-1 to 2:
			for ix = nx-1:-1:2
				T(ix) = (d(ix) + A1(M(ix),M(ix+1))*T(ix+1))/b(ix);
			end
			% Back sweep, index 1:
			if ~heatsinkedboundary
				T(1) = (d(1) + A1(M(1),M(2))*T(2))/b(1);
			end
			
			h_line.YData = T;
            h_title.String = num2str(it);
            drawnow;
        end
    case 2
		%% Initialize arrays
		x = dx*(-(nx-1)/2:(nx-1)/2);
		y = dy*(-(ny-1)/2:(ny-1)/2);
        N = nx*ny;
        M = ones(nx,ny);
        M(5:15,23:30) = 2;
        ABS = zeros(nx,ny); % Power deposited per unit length/area/volume per unit time
        ABS(5,20) = 100;
		ABS(17,25) = 50;
        
		A1 = NaN(nM,nM);
        A2 = NaN(nM,nM);
		s = zeros(nx,ny);
		for i=1:nM
			for j=1:nM
                TC_eff = 2*TC(i)*TC(j)/(TC(i)+TC(j));
                if isnan(TC_eff); TC_eff = 0; end
                A1(i,j) = dt*TC_eff/(2*VHC(i)*dx^2); % A1 coupling factor for heat flow along x from medium j into medium i
                A2(i,j) = dt*TC_eff/(2*VHC(i)*dy^2); % A2 coupling factor for heat flow along y from medium j into medium i
			end
		end
		for iy = 1:ny
			for ix = 1:nx
				heatsinkedvoxel = heatsinkedboundary && (ix == 1 || ix == nx || iy == 1 || iy == ny);
				if ~heatsinkedvoxel
					s(ix,iy) = ABS(ix,iy)*dt/VHC(M(ix,iy));
				end
			end
		end
        
		%% Initialize figures
        figure(2);clf;
        imagesc(x,y,M.');
		colorbar;
        axis equal
        axis tight
        axis xy
        
        T = zeros(nx,ny);
        d = zeros(nx,ny);
		b = zeros(max(nx,ny),1);
        figure(1);clf;
        h_im = imagesc(x,y,T.');
        colorbar;
        axis equal
        axis tight
        axis xy
        h_title = title('0');

		%% Do loop
        for it=1:nt
			%% Explicit part of step 1
			for iy = 1:ny
				for ix = 1:nx % construct d, the right hand side of the system of linear equations
					if ~(heatsinkedboundary && (ix == 1 || ix == nx || iy == 1 || iy == ny))
						delta = 0;
						if ix ~= 1;  delta = delta +   A1(M(ix,iy),M(ix-1,iy))*(T(ix-1,iy)-T(ix,iy)); end
						if ix ~= nx; delta = delta +   A1(M(ix,iy),M(ix+1,iy))*(T(ix+1,iy)-T(ix,iy)); end
						if iy ~= 1;  delta = delta + 2*A2(M(ix,iy),M(ix,iy-1))*(T(ix,iy-1)-T(ix,iy)); end
						if iy ~= ny; delta = delta + 2*A2(M(ix,iy),M(ix,iy+1))*(T(ix,iy+1)-T(ix,iy)); end
					else
						delta = 0;
					end
					d(ix,iy) = T(ix,iy) + delta;
				end
			end
			
			%% Implicit part of step 1, Thomas algorithm along x, sweeps up from 2 to nx and then down from nx to 1
			for iy = (1 + heatsinkedboundary):(ny - heatsinkedboundary) % heatsinkedboundary may be 0 or 1
				% Forward sweep, index 2:
				if heatsinkedboundary
					b(1) = 1;
				else
					b(1) = 1 + A1(M(1,iy),M(2,iy));
				end
				w = -A1(M(2,iy),M(1,iy))/b(1);
				if heatsinkedboundary
					b(2) = 1 + A1(M(2,iy),M(3,iy)) + A1(M(2,iy),M(1,iy));
				else
					b(2) = 1 + A1(M(2,iy),M(3,iy)) + A1(M(2,iy),M(1,iy)) + w*A1(M(1,iy),M(2,iy));
				end
				d(2,iy) = d(2,iy) - w*d(1,iy);
				% Forward sweep, indices 3 to nx-1:
				for ix = 3:nx-1
					w = -A1(M(ix,iy),M(ix-1,iy))/b(ix-1);
					b(ix) = 1 + A1(M(ix,iy),M(ix-1,iy)) + A1(M(ix,iy),M(ix+1,iy));
					b(ix) = b(ix) + w*A1(M(ix-1,iy),M(ix,iy));
					d(ix,iy) = d(ix,iy) - w*d(ix-1,iy);
				end
				% Forward sweep, index nx:
				if heatsinkedboundary
					b(nx) = 1;
				else
					w = -A1(M(nx,iy),M(nx-1,iy))/b(nx-1);
					b(nx) = 1 + A1(M(nx,iy),M(nx-1,iy));
					b(nx) = b(nx) + w*A1(M(nx-1,iy),M(nx,iy));
					d(nx,iy) = d(nx,iy) - w*d(nx-1,iy);
				end

				% Back sweep, index nx:
				d(nx,iy) = d(nx,iy)/b(nx);
				% Back sweep, indices nx-1 to 2:
				for ix = nx-1:-1:2
					d(ix,iy) = (d(ix,iy) + A1(M(ix,iy),M(ix+1,iy))*d(ix+1,iy))/b(ix);
				end
				% Back sweep, index 1:
				if ~heatsinkedboundary
					d(1,iy) = (d(1,iy) + A1(M(1,iy),M(2,iy))*d(2,iy))/b(1);
				end
			end
			
			%% Explicit part of step 2
			for iy = 1:ny
				for ix = 1:nx % construct d, the right hand side of the system of linear equations
					if ~(heatsinkedboundary && (ix == 1 || ix == nx || iy == 1 || iy == ny))
						delta = 0;
						if iy ~= 1;  delta = delta - A2(M(ix,iy),M(ix,iy-1))*(T(ix,iy-1)-T(ix,iy)); end
						if iy ~= ny; delta = delta - A2(M(ix,iy),M(ix,iy+1))*(T(ix,iy+1)-T(ix,iy)); end
					else
						delta = 0;
					end
					d(ix,iy) = d(ix,iy) + delta;
				end
			end
			
			%% Implicit part of step 2, Thomas algorithm along y, sweeps up from 2 to ny and then down from ny to 1
			for ix = (1 + heatsinkedboundary):(nx - heatsinkedboundary) % heatsinkedboundary may be 0 or 1
				% Forward sweep, index 2:
				if heatsinkedboundary
					b(1) = 1;
				else
					b(1) = 1 + A2(M(ix,1),M(ix,2));
				end
				w = -A2(M(ix,2),M(ix,1))/b(1);
				if heatsinkedboundary
					b(2) = 1 + A2(M(ix,2),M(ix,3)) + A2(M(ix,2),M(ix,1));
				else
					b(2) = 1 + A2(M(ix,2),M(ix,3)) + A2(M(ix,2),M(ix,1)) + w*A2(M(ix,1),M(ix,2));
				end
				d(ix,2) = d(ix,2) - w*d(ix,1);
				% Forward sweep, indices 3 to ny-1:
				for iy = 3:ny-1
					w = -A2(M(ix,iy),M(ix,iy-1))/b(iy-1);
					b(iy) = 1 + A2(M(ix,iy),M(ix,iy-1)) + A2(M(ix,iy),M(ix,iy+1));
					b(iy) = b(iy) + w*A2(M(ix,iy-1),M(ix,iy));
					d(ix,iy) = d(ix,iy) - w*d(ix,iy-1);
				end
				% Forward sweep, index ny:
				if heatsinkedboundary
					b(ny) = 1;
				else
					w = -A2(M(ix,ny),M(ix,ny-1))/b(ny-1);
					b(ny) = 1 + A2(M(ix,ny),M(ix,ny-1));
					b(ny) = b(ny) + w*A2(M(ix,ny-1),M(ix,ny));
					d(ix,ny) = d(ix,ny) - w*d(ix,ny-1);
				end

				% Back sweep, index ny:
				T(ix,ny) = d(ix,ny)/b(ny);
				% Back sweep, indices ny-1 to 2:
				for iy = ny-1:-1:2
					T(ix,iy) = (d(ix,iy) + A2(M(ix,iy),M(ix,iy+1))*T(ix,iy+1))/b(iy);
					T(ix,iy+1) = T(ix,iy+1) + s(ix,iy+1);
				end
				% Back sweep, index 1:
				if ~heatsinkedboundary
					T(ix,1) = (d(ix,1) + A2(M(ix,1),M(ix,2))*T(ix,2))/b(1);
				end
				T(ix,2) = T(ix,2) + s(ix,2);
				T(ix,1) = T(ix,1) + s(ix,1);
			end
			
            h_im.CData = reshape(T,[nx ny]).';
            h_title.String = num2str(it);
            drawnow;
		end
	case 3
		%% Initialize arrays
		x = dx*(-(nx-1)/2:(nx-1)/2);
		y = dy*(-(ny-1)/2:(ny-1)/2);
		z = dz*(0.5:nz-0.5);
        N = nx*ny*nz;
        M = ones(nx,ny,nz);
        M(1:end-6,:,10:end) = 2;
        ABS = zeros(nx,ny,nz); % Power deposited per unit length/area/volume per unit time
%         ABS(5,20,5) = 100;
% 		ABS(17,25,10) = 50;
% 		ABS(round((nx-1)/2+1),round((ny-1)/2+1),round((nz-1)/2+1)) = 100;
        
		A1 = NaN(nM,nM);
        A2 = NaN(nM,nM);
		A3 = NaN(nM,nM);
		s = zeros(nx,ny,nz);
		for i=1:nM
			for j=1:nM
                TC_eff = 2*TC(i)*TC(j)/(TC(i)+TC(j));
                if isnan(TC_eff); TC_eff = 0; end
                A1(i,j) = dt*TC_eff/(2*VHC(i)*dx^2); % A1 coupling factor for heat flow along x from medium j into medium i
                A2(i,j) = dt*TC_eff/(2*VHC(i)*dy^2); % A2 coupling factor for heat flow along y from medium j into medium i
                A3(i,j) = dt*TC_eff/(2*VHC(i)*dz^2); % A3 coupling factor for heat flow along z from medium j into medium i
			end
		end
		for iz = 1:nz
			for iy = 1:ny
				for ix = 1:nx
					heatsinkedvoxel = heatsinkedboundary && (ix == 1 || ix == nx || iy == 1 || iy == ny || iz == 1 || iz == nz);
					if ~heatsinkedvoxel
						s(ix,iy,iz) = ABS(ix,iy,iz)*dt/VHC(M(ix,iy,iz));
					end
				end
			end
		end
		
		%% Initialize figures
		plotVolumetric(2,x,y,z,M,'MCmatlab_heat');
        axis equal
        axis tight
        
        T = zeros(nx,ny,nz);
		T(16,16,16) = 1;
        d = zeros(nx,ny,nz);
		b = zeros(max([nx ny nz]),1);
        h_heatfig = plotVolumetric(1,x,y,z,T,'MCmatlab_heat','slicePositions',[0.5 0.5 0.5]);
        h_title = title('0');
		caxis([0 0.02]);

		%% Do loop
		for it=1:nt
			%% Explicit part of step 1
			for iz = 1:nz
				for iy = 1:ny
					for ix = 1:nx % construct d, the right hand side of the system of linear equations
						if ~(heatsinkedboundary && (ix == 1 || ix == nx || iy == 1 || iy == ny || iz == 1 || iz == nz))
							delta = 0;
							if ix ~= 1;  delta = delta +   A1(M(ix,iy,iz),M(ix-1,iy,iz))*(T(ix-1,iy,iz)-T(ix,iy,iz)); end
							if ix ~= nx; delta = delta +   A1(M(ix,iy,iz),M(ix+1,iy,iz))*(T(ix+1,iy,iz)-T(ix,iy,iz)); end
							if iy ~= 1;  delta = delta + 2*A2(M(ix,iy,iz),M(ix,iy-1,iz))*(T(ix,iy-1,iz)-T(ix,iy,iz)); end
							if iy ~= ny; delta = delta + 2*A2(M(ix,iy,iz),M(ix,iy+1,iz))*(T(ix,iy+1,iz)-T(ix,iy,iz)); end
							if iz ~= 1;  delta = delta + 2*A3(M(ix,iy,iz),M(ix,iy,iz-1))*(T(ix,iy,iz-1)-T(ix,iy,iz)); end
							if iz ~= nz; delta = delta + 2*A3(M(ix,iy,iz),M(ix,iy,iz+1))*(T(ix,iy,iz+1)-T(ix,iy,iz)); end
						else
							delta = 0;
						end
						d(ix,iy,iz) = T(ix,iy,iz) + delta;
					end
				end
			end
			
			%% Implicit part of step 1, Thomas algorithm along x, sweeps up from 2 to nx and then down from nx to 1
			for iz = (1 + heatsinkedboundary):(nz - heatsinkedboundary) % heatsinkedboundary may be 0 or 1
				for iy = (1 + heatsinkedboundary):(ny - heatsinkedboundary) % heatsinkedboundary may be 0 or 1
					% Forward sweep, index 2:
					if heatsinkedboundary
						b(1) = 1;
					else
						b(1) = 1 + A1(M(1,iy,iz),M(2,iy,iz));
					end
					w = -A1(M(2,iy,iz),M(1,iy,iz))/b(1);
					if heatsinkedboundary
						b(2) = 1 + A1(M(2,iy,iz),M(3,iy,iz)) + A1(M(2,iy,iz),M(1,iy,iz));
					else
						b(2) = 1 + A1(M(2,iy,iz),M(3,iy,iz)) + A1(M(2,iy,iz),M(1,iy,iz)) + w*A1(M(1,iy,iz),M(2,iy,iz));
					end
					d(2,iy,iz) = d(2,iy,iz) - w*d(1,iy,iz);
					% Forward sweep, indices 3 to nx-1:
					for ix = 3:nx-1
						w = -A1(M(ix,iy,iz),M(ix-1,iy,iz))/b(ix-1);
						b(ix) = 1 + A1(M(ix,iy,iz),M(ix-1,iy,iz)) + A1(M(ix,iy,iz),M(ix+1,iy,iz));
						b(ix) = b(ix) + w*A1(M(ix-1,iy,iz),M(ix,iy,iz));
						d(ix,iy,iz) = d(ix,iy,iz) - w*d(ix-1,iy,iz);
					end
					% Forward sweep, index nx:
					if heatsinkedboundary
						b(nx) = 1;
					else
						w = -A1(M(nx,iy,iz),M(nx-1,iy,iz))/b(nx-1);
						b(nx) = 1 + A1(M(nx,iy,iz),M(nx-1,iy,iz));
						b(nx) = b(nx) + w*A1(M(nx-1,iy,iz),M(nx,iy,iz));
						d(nx,iy,iz) = d(nx,iy,iz) - w*d(nx-1,iy,iz);
					end

					% Back sweep, index nx:
					d(nx,iy,iz) = d(nx,iy,iz)/b(nx);
					% Back sweep, indices nx-1 to 2:
					for ix = nx-1:-1:2
						d(ix,iy,iz) = (d(ix,iy,iz) + A1(M(ix,iy,iz),M(ix+1,iy,iz))*d(ix+1,iy,iz))/b(ix);
					end
					% Back sweep, index 1:
					if ~heatsinkedboundary
						d(1,iy,iz) = (d(1,iy,iz) + A1(M(1,iy,iz),M(2,iy,iz))*d(2,iy,iz))/b(1);
					end
				end
			end
			
			%% Explicit part of step 2
			for iz = 1:nz
				for iy = 1:ny
					for ix = 1:nx % construct d, the right hand side of the system of linear equations
						if ~(heatsinkedboundary && (ix == 1 || ix == nx || iy == 1 || iy == ny || iz == 1 || iz == nz))
							delta = 0;
							if iy ~= 1;  delta = delta - A2(M(ix,iy,iz),M(ix,iy-1,iz))*(T(ix,iy-1,iz)-T(ix,iy,iz)); end
							if iy ~= ny; delta = delta - A2(M(ix,iy,iz),M(ix,iy+1,iz))*(T(ix,iy+1,iz)-T(ix,iy,iz)); end
						else
							delta = 0;
						end
						d(ix,iy,iz) = d(ix,iy,iz) + delta;
					end
				end
			end
			
			%% Implicit part of step 2, Thomas algorithm along y, sweeps up from 2 to ny and then down from ny to 1
			for iz = (1 + heatsinkedboundary):(nz - heatsinkedboundary) % heatsinkedboundary may be 0 or 1
				for ix = (1 + heatsinkedboundary):(nx - heatsinkedboundary) % heatsinkedboundary may be 0 or 1
					% Forward sweep, index 2:
					if heatsinkedboundary
						b(1) = 1;
					else
						b(1) = 1 + A2(M(ix,1,iz),M(ix,2,iz));
					end
					w = -A2(M(ix,2,iz),M(ix,1,iz))/b(1);
					if heatsinkedboundary
						b(2) = 1 + A2(M(ix,2,iz),M(ix,3,iz)) + A2(M(ix,2,iz),M(ix,1,iz));
					else
						b(2) = 1 + A2(M(ix,2,iz),M(ix,3,iz)) + A2(M(ix,2,iz),M(ix,1,iz)) + w*A2(M(ix,1,iz),M(ix,2,iz));
					end
					d(ix,2,iz) = d(ix,2,iz) - w*d(ix,1,iz);
					% Forward sweep, indices 3 to ny-1:
					for iy = 3:ny-1
						w = -A2(M(ix,iy,iz),M(ix,iy-1,iz))/b(iy-1);
						b(iy) = 1 + A2(M(ix,iy,iz),M(ix,iy-1,iz)) + A2(M(ix,iy,iz),M(ix,iy+1,iz));
						b(iy) = b(iy) + w*A2(M(ix,iy-1,iz),M(ix,iy,iz));
						d(ix,iy,iz) = d(ix,iy,iz) - w*d(ix,iy-1,iz);
					end
					% Forward sweep, index ny:
					if heatsinkedboundary
						b(ny) = 1;
					else
						w = -A2(M(ix,ny,iz),M(ix,ny-1,iz))/b(ny-1);
						b(ny) = 1 + A2(M(ix,ny,iz),M(ix,ny-1,iz));
						b(ny) = b(ny) + w*A2(M(ix,ny-1,iz),M(ix,ny,iz));
						d(ix,ny,iz) = d(ix,ny,iz) - w*d(ix,ny-1,iz);
					end

					% Back sweep, index ny:
					d(ix,ny,iz) = d(ix,ny,iz)/b(ny);
					% Back sweep, indices ny-1 to 2:
					for iy = ny-1:-1:2
						d(ix,iy,iz) = (d(ix,iy,iz) + A2(M(ix,iy,iz),M(ix,iy+1,iz))*d(ix,iy+1,iz))/b(iy);
					end
					% Back sweep, index 1:
					if ~heatsinkedboundary
						d(ix,1,iz) = (d(ix,1,iz) + A2(M(ix,1,iz),M(ix,2,iz))*d(ix,2,iz))/b(1);
					end
				end
			end
			
			%% Explicit part of step 3
			for iz = 1:nz
				for iy = 1:ny
					for ix = 1:nx % construct d, the right hand side of the system of linear equations
						if ~(heatsinkedboundary && (ix == 1 || ix == nx || iy == 1 || iy == ny || iz == 1 || iz == nz))
							delta = 0;
							if iz ~= 1;  delta = delta - A3(M(ix,iy,iz),M(ix,iy,iz-1))*(T(ix,iy,iz-1)-T(ix,iy,iz)); end
							if iz ~= nz; delta = delta - A3(M(ix,iy,iz),M(ix,iy,iz+1))*(T(ix,iy,iz+1)-T(ix,iy,iz)); end
						else
							delta = 0;
						end
						d(ix,iy,iz) = d(ix,iy,iz) + delta;
					end
				end
			end
			
			%% Implicit part of step 3, Thomas algorithm along z, sweeps up from 2 to nz and then down from nz to 1
			for iy = (1 + heatsinkedboundary):(ny - heatsinkedboundary) % heatsinkedboundary may be 0 or 1
				for ix = (1 + heatsinkedboundary):(nx - heatsinkedboundary) % heatsinkedboundary may be 0 or 1
					% Forward sweep, index 2:
					if heatsinkedboundary
						b(1) = 1;
					else
						b(1) = 1 + A3(M(ix,iy,1),M(ix,iy,2));
					end
					w = -A3(M(ix,iy,2),M(ix,iy,1))/b(1);
					if heatsinkedboundary
						b(2) = 1 + A3(M(ix,iy,2),M(ix,iy,3)) + A3(M(ix,iy,2),M(ix,iy,1));
					else
						b(2) = 1 + A3(M(ix,iy,2),M(ix,iy,3)) + A3(M(ix,iy,2),M(ix,iy,1)) + w*A3(M(ix,iy,1),M(ix,iy,2));
					end
					d(ix,iy,2) = d(ix,iy,2) - w*d(ix,iy,1);
					% Forward sweep, indices 3 to nz-1:
					for iz = 3:nz-1
						w = -A3(M(ix,iy,iz),M(ix,iy,iz-1))/b(iz-1);
						b(iz) = 1 + A3(M(ix,iy,iz),M(ix,iy,iz-1)) + A3(M(ix,iy,iz),M(ix,iy,iz+1));
						b(iz) = b(iz) + w*A3(M(ix,iy,iz-1),M(ix,iy,iz));
						d(ix,iy,iz) = d(ix,iy,iz) - w*d(ix,iy,iz-1);
					end
					% Forward sweep, index nz:
					if heatsinkedboundary
						b(nz) = 1;
					else
						w = -A3(M(ix,iy,nz),M(ix,iy,nz-1))/b(nz-1);
						b(nz) = 1 + A3(M(ix,iy,nz),M(ix,iy,nz-1));
						b(nz) = b(nz) + w*A3(M(ix,iy,nz-1),M(ix,iy,nz));
						d(ix,iy,nz) = d(ix,iy,nz) - w*d(ix,iy,nz-1);
					end

					% Back sweep, index nz:
					T(ix,iy,nz) = d(ix,iy,nz)/b(nz);
					% Back sweep, indices nz-1 to 2:
					for iz = nz-1:-1:2
						T(ix,iy,iz) = (d(ix,iy,iz) + A3(M(ix,iy,iz),M(ix,iy,iz+1))*T(ix,iy,iz+1))/b(iz);
						T(ix,iy,iz+1) = T(ix,iy,iz+1) + s(ix,iy,iz+1); % Heat deposition (not strictly speaking a part of the implicit step)
					end
					% Back sweep, index 1:
					if ~heatsinkedboundary
						T(ix,iy,1) = (d(ix,iy,1) + A3(M(ix,iy,1),M(ix,iy,2))*T(ix,iy,2))/b(1);
					end
					T(ix,iy,2) = T(ix,iy,2) + s(ix,iy,2); % Heat deposition (not strictly speaking a part of the implicit step)
					T(ix,iy,1) = T(ix,iy,1) + s(ix,iy,1); % Heat deposition (not strictly speaking a part of the implicit step)
				end
			end
			
			updateVolumetric(h_heatfig,T);
            h_title.String = num2str(it);
            drawnow;
			if it == 1
				pause
			else
				pause(0.1)
			end
		end
	case 4
		clear Ginput
		Ginput.matchedInterfaces = true; % Assumes all refractive indices are 1
		Ginput.boundaryType      = 1; % 0: No escaping boundaries, 1: All cuboid boundaries are escaping, 2: Top cuboid boundary only is escaping

		Ginput.wavelength        = 532; % [nm] Excitation wavelength, used for determination of optical properties for excitation light

		Ginput.nx                = nx; % Number of bins in the x direction
		Ginput.ny                = ny; % Number of bins in the y direction
		Ginput.nz                = nz; % Number of bins in the z direction
		Ginput.Lx                = nx*dx; % [cm] x size of simulation cuboid
		Ginput.Ly                = ny*dy; % [cm] y size of simulation cuboid
		Ginput.Lz                = nz*dz; % [cm] z size of simulation cuboid

		Ginput.GeomFunc          = @GeometryDefinition_StandardTissue; % Function to use for defining the distribution of media in the cuboid. Defined at the end of this m file.

		% Execution, do not modify the next two lines:
		G = defineGeometry(Ginput);
		
		x = G.x;
		y = G.y;
		z = G.z;
		nM = length(G.mediaProperties); % Number of different media in simulation

		dTdt_abs = zeros(G.nx,G.ny,G.nz); % Power deposited per unit length/area/volume per unit time
%         ABS(5,20,5) = 100;
% 		ABS(17,25,10) = 50;
% 		dTdt_abs(round((G.nx-1)/2+1),round((G.ny-1)/2+1),round((G.nz-1)/2+1)) = 100;

		VHC = [G.mediaProperties.VHC]; % Volumetric heat capacity array. Element i is the VHC of medium i in the reduced mediaProperties list.
		TC = [G.mediaProperties.TC]; % Thermal conductivity array. Element i is the TC of medium i in the reduced mediaProperties list.
		A = [G.mediaProperties.A];
		E = [G.mediaProperties.E];
		if(any(A)) % If non-zero Arrhenius data exists, prepare to calculate thermal damage.
			HSoutput.Omega = zeros(size(G.M));
		else
			HSoutput.Omega = NaN;
		end
		TC_eff = 2*(TC'*TC)./(ones(length(TC))*diag(TC)+diag(TC)*ones(length(TC))); % Same as TC_eff(i,j) = 2*TC_red(i)*TC_red(j)/(TC_red(i)+TC_red(j)) but without for loops
		TC_eff(isnan(TC_eff)) = 0; % Neighboring insulating voxels return NaN but should just be 0

		dTdtperdeltaT  = cat(3,TC_eff./(diag(VHC)*ones(nM))/G.dx^2,...
							   TC_eff./(diag(VHC)*ones(nM))/G.dy^2,...
							   TC_eff./(diag(VHC)*ones(nM))/G.dz^2); % [deg C /(s*deg C)] Time derivative of voxel temperature per voxel-to-voxel temperature difference. Third dimension corresponds to the different directions of heat diffusion (x, y and z)
		clear TC_eff

		heatSimParameters = struct('M',G.M-1,'A',A,'E',E,'dTdtperdeltaT',dTdtperdeltaT,'dTdt_abs',dTdt_abs,...
                'useAllCPUs',true,'heatBoundaryType',double(heatsinkedboundary),...
                'tempSensorCornerIdxs',[],'tempSensorInterpWeights',[]); % Contents of G.M have to be converted from Matlab's 1-based indexing to C's 0-based indexing.
        heatSimParameters.lightsOn = true;
        heatSimParameters.steps = 1;
		heatSimParameters.dt = dt;

        plotVolumetric(1,G.x,G.y,G.z,G.M,'MCmatlab_GeometryIllustration',G.mediaProperties,'slicePositions',[1 1 0.5]);
        
		T = ones(G.nx,G.ny,G.nz);
% 		T(16,16,16) = 1;

        h_heatfig = plotVolumetric(2,G.x,G.y,G.z,T,'MCmatlab_heat','slicePositions',[0.5 0.5 0.5]);
        h_title = title('1');
		caxis([0 1]);
		clear totalenergy
		for it=1:nt
			[T,HSoutput.Omega,newSensorTemps] = finiteElementHeatPropagator(T,HSoutput.Omega,heatSimParameters);
			updateVolumetric(h_heatfig,T);
			totalenergy(it) = sum(VHC(G.M(:)).*T(:)*G.dx*G.dy*G.dz);
			h_title.String = num2str(it);
		end
		figure;plot(totalenergy);
end

function M = GeometryDefinition_StandardTissue(X,Y,Z,parameters)
tissuedepth = 0.03;
M = ones(size(X)); % Air
M(Z > tissuedepth) = 3; % "Standard" tissue
M(X > 0.03) = 1;
M(:) = 3;
end