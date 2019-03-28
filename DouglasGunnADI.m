%% User specified parameters
calctype = 2; % 1: 1D, 2:2D, 3:3D

heatsinkedboundary = false;

dt = 0.0001;

dx = 0.1;
dy = 0.1;
dz = 0.08;

TC = [10 ; 1];

VHC = [1 ; 2];

nx = 40;
ny = 40;
nz = 30;
nt = 500;

%% Calculations
nM = length(TC);
x = dx*(-(nx-1)/2:(nx-1)/2);
y = dy*(-(ny-1)/2:(ny-1)/2);
switch calctype
    case 1
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
			xidx = i;
            heatsinkedvoxel = heatsinkedboundary && (xidx == 1 || xidx == nx);
			if ~heatsinkedvoxel
				s(i) = ABS(i)*dt/VHC(M(i));
			end
		end

        T = zeros(N,1);
		d = zeros(N,1);
		T(20) = 1;
		b = NaN(nx,1);
		figure(1);clf;
        h_line = plot(T);
        h_title = title('0');
        for tidx=1:nt
			for xidx = 1:nx % construct d, the right hand side of the system of linear equations
				if ~(heatsinkedboundary && (xidx == 1 || xidx == nx))
					delta = s(xidx);
					if xidx ~= 1;  delta = delta + A1(M(xidx),M(xidx-1))*(T(xidx-1)-T(xidx)); end
					if xidx ~= nx; delta = delta + A1(M(xidx),M(xidx+1))*(T(xidx+1)-T(xidx)); end
				else
					delta = 0;
				end
				d(xidx) = T(xidx) + delta;
			end
			
			%% Thomas algorithm, sweeps up from 2 to nx and then down from nx to 1
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
			for xidx = 3:nx-1
				w = -A1(M(xidx),M(xidx-1))/b(xidx-1);
				b(xidx) = 1 + A1(M(xidx),M(xidx-1)) + A1(M(xidx),M(xidx+1));
				b(xidx) = b(xidx) + w*A1(M(xidx-1),M(xidx));
				d(xidx) = d(xidx) - w*d(xidx-1);
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
			for xidx = nx-1:-1:2
				T(xidx) = (d(xidx) + A1(M(xidx),M(xidx+1))*T(xidx+1))/b(xidx);
			end
			% Back sweep, index 1:
			if ~heatsinkedboundary
				T(1) = (d(1) + A1(M(1),M(2))*T(2))/b(1);
			end
			
			h_line.YData = T;
            h_title.String = num2str(tidx);
            drawnow;
        end
        
        
        
        
        
        
        
        
        
        
    case 2
        N = nx*ny;
        M = ones(nx,ny);
        M(5:15,23:30) = 2;
        ABS = zeros(nx,ny); % Power deposited per unit length/area/volume per unit time
        ABS(5,20) = 100;
        
		A1 = NaN(nM,nM);
        A2 = NaN(nM,nM);
		s = zeros(nx,ny);
		for i=1:nM
			for j=1:nM
                TC_eff = 2*TC(i)*TC(j)/(TC(i)+TC(j));
                if isnan(TC_eff); TC_eff = 0; end
                A1(i,j) = dt/2*TC_eff/(2*VHC(i)*dx^2); % A1 coupling factor for heat flow along x from medium j into medium i
                A2(i,j) = dt/2*TC_eff/(2*VHC(i)*dy^2); % A2 coupling factor for heat flow along y from medium j into medium i
			end
		end
		for xidx = 1:nx
			for yidx = 1:ny
				heatsinkedvoxel = heatsinkedboundary && (xidx == 1 || xidx == nx || yidx == 1 || yidx == ny);
				if ~heatsinkedvoxel
					s(xidx,yidx) = ABS(xidx,yidx)*dt/VHC(M(xidx,yidx));
				end
			end
		end
        
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
        h_im = imagesc(x,y,reshape(T,[nx ny]).');
        colorbar;
        axis equal
        axis tight
        axis xy
        h_title = title('0');

        for tidx=1:nt
			%% Explicit part of step 1
			for yidx = 1:ny
				for xidx = 1:nx % construct d, the right hand side of the system of linear equations
					if ~(heatsinkedboundary && (xidx == 1 || xidx == nx || yidx == 1 || yidx == ny))
						delta = s(xidx,yidx);
						if xidx ~= 1;  delta = delta +   A1(M(xidx,yidx),M(xidx-1,yidx))*(T(xidx-1,yidx)-T(xidx,yidx)); end
						if xidx ~= nx; delta = delta +   A1(M(xidx,yidx),M(xidx+1,yidx))*(T(xidx+1,yidx)-T(xidx,yidx)); end
						if yidx ~= 1;  delta = delta + 2*A2(M(xidx,yidx),M(xidx,yidx-1))*(T(xidx,yidx-1)-T(xidx,yidx)); end
						if yidx ~= ny; delta = delta + 2*A2(M(xidx,yidx),M(xidx,yidx+1))*(T(xidx,yidx+1)-T(xidx,yidx)); end
					else
						delta = 0;
					end
					d(xidx,yidx) = T(xidx,yidx) + delta;
				end
			end
			
			%% Implicit part of step 1, Thomas algorithm along x, sweeps up from 2 to nx and then down from nx to 1
			for yidx = (1 + heatsinkedboundary):(ny - heatsinkedboundary) % heatsinkedboundary may be 0 or 1
				% Forward sweep, index 2:
				if heatsinkedboundary
					b(1) = 1;
				else
					b(1) = 1 + A1(M(1,yidx),M(2,yidx));
				end
				w = -A1(M(2,yidx),M(1,yidx))/b(1);
				if heatsinkedboundary
					b(2) = 1 + A1(M(2,yidx),M(3,yidx)) + A1(M(2,yidx),M(1,yidx));
				else
					b(2) = 1 + A1(M(2,yidx),M(3,yidx)) + A1(M(2,yidx),M(1,yidx)) + w*A1(M(1,yidx),M(2,yidx));
				end
				d(2,yidx) = d(2,yidx) - w*d(1,yidx);
				% Forward sweep, indices 3 to nx-1:
				for xidx = 3:nx-1
					w = -A1(M(xidx,yidx),M(xidx-1,yidx))/b(xidx-1);
					b(xidx) = 1 + A1(M(xidx,yidx),M(xidx-1,yidx)) + A1(M(xidx,yidx),M(xidx+1,yidx));
					b(xidx) = b(xidx) + w*A1(M(xidx-1,yidx),M(xidx,yidx));
					d(xidx,yidx) = d(xidx,yidx) - w*d(xidx-1,yidx);
				end
				% Forward sweep, index nx:
				if heatsinkedboundary
					b(nx) = 1;
				else
					w = -A1(M(nx,yidx),M(nx-1,yidx))/b(nx-1);
					b(nx) = 1 + A1(M(nx,yidx),M(nx-1,yidx));
					b(nx) = b(nx) + w*A1(M(nx-1,yidx),M(nx,yidx));
					d(nx,yidx) = d(nx,yidx) - w*d(nx-1,yidx);
				end

				% Back sweep, index nx:
				d(nx,yidx) = d(nx,yidx)/b(nx);
				% Back sweep, indices nx-1 to 2:
				for xidx = nx-1:-1:2
					d(xidx,yidx) = (d(xidx,yidx) + A1(M(xidx,yidx),M(xidx+1,yidx))*d(xidx+1,yidx))/b(xidx);
				end
				% Back sweep, index 1:
				if ~heatsinkedboundary
					d(1,yidx) = (d(1,yidx) + A1(M(1,yidx),M(2,yidx))*d(2,yidx))/b(1);
				end
			end
			
			%% Explicit part of step 2
			for yidx = 1:ny
				for xidx = 1:nx % construct d, the right hand side of the system of linear equations
					if ~(heatsinkedboundary && (xidx == 1 || xidx == nx || yidx == 1 || yidx == ny))
						delta = 0;
						if yidx ~= 1;  delta = delta - A2(M(xidx,yidx),M(xidx,yidx-1))*(T(xidx,yidx-1)-T(xidx,yidx)); end
						if yidx ~= ny; delta = delta - A2(M(xidx,yidx),M(xidx,yidx+1))*(T(xidx,yidx+1)-T(xidx,yidx)); end
					else
						delta = 0;
					end
					d(xidx,yidx) = d(xidx,yidx) + delta;
				end
			end
			
			%% Implicit part of step 2, Thomas algorithm along y, sweeps up from 2 to ny and then down from ny to 1
			for xidx = (1 + heatsinkedboundary):(nx - heatsinkedboundary) % heatsinkedboundary may be 0 or 1
				% Forward sweep, index 2:
				if heatsinkedboundary
					b(1) = 1;
				else
					b(1) = 1 + A2(M(xidx,1),M(xidx,2));
				end
				w = -A2(M(xidx,2),M(xidx,1))/b(1);
				if heatsinkedboundary
					b(2) = 1 + A2(M(xidx,2),M(xidx,3)) + A2(M(xidx,2),M(xidx,1));
				else
					b(2) = 1 + A2(M(xidx,2),M(xidx,3)) + A2(M(xidx,2),M(xidx,1)) + w*A2(M(xidx,1),M(xidx,2));
				end
				d(xidx,2) = d(xidx,2) - w*d(xidx,1);
				% Forward sweep, indices 3 to ny-1:
				for yidx = 3:ny-1
					w = -A2(M(xidx,yidx),M(xidx,yidx-1))/b(yidx-1);
					b(yidx) = 1 + A2(M(xidx,yidx),M(xidx,yidx-1)) + A2(M(xidx,yidx),M(xidx,yidx+1));
					b(yidx) = b(yidx) + w*A2(M(xidx,yidx-1),M(xidx,yidx));
					d(xidx,yidx) = d(xidx,yidx) - w*d(xidx,yidx-1);
				end
				% Forward sweep, index ny:
				if heatsinkedboundary
					b(ny) = 1;
				else
					w = -A2(M(xidx,ny),M(xidx,ny-1))/b(ny-1);
					b(ny) = 1 + A2(M(xidx,ny),M(xidx,ny-1));
					b(ny) = b(ny) + w*A2(M(xidx,ny-1),M(xidx,ny));
					d(xidx,ny) = d(xidx,ny) - w*d(xidx,ny-1);
				end

				% Back sweep, index ny:
				T(xidx,ny) = d(xidx,ny)/b(ny);
				% Back sweep, indices ny-1 to 2:
				for yidx = ny-1:-1:2
					T(xidx,yidx) = (d(xidx,yidx) + A2(M(xidx,yidx),M(xidx,yidx+1))*T(xidx,yidx+1))/b(yidx);
				end
				% Back sweep, index 1:
				if ~heatsinkedboundary
					T(xidx,1) = (d(xidx,1) + A2(M(xidx,1),M(xidx,2))*T(xidx,2))/b(1);
				end
			end
			
            h_im.CData = reshape(T,[nx ny]).';
            h_title.String = num2str(tidx);
            drawnow;
        end
end
