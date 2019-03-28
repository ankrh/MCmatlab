calctype = 1; % 1: 1D, 2:2D, 3:3D

heatsinkedboundary = false;

dt = 0.0001;

dx = 0.1;
dy = 0.03;
dz = 0.08;

TC = [10 ; 1];

VHC = [1 ; 2];

nx = 20;
ny = 10;
nz = 30;
nt = 500;

%% Calculate matrices used for Douglas-Gunn ADI method
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
		B = NaN(nM,nM,1);
		s = zeros(nx,1);
		for i=1:nM
			for j=1:nM
				TC_eff = 2*TC(i)*TC(j)/(TC(i)+TC(j));
				A1(i,j) = dt*TC_eff/(2*VHC(i)*dx^2); % A1 coupling factor for heat flow from medium j into medium i
				B(i,j)  = dt*TC_eff/(2*VHC(i)*dx^2); % B coupling factor for heat flow from medium j into medium i
			end
		end
		for i=1:N
            heatsinkedvoxel = heatsinkedboundary && (ix == 1 || ix == nx);
            if ~heatsinkedvoxel
                s(i) = ABS(i)*dt/VHC(M(i));
			end
		end

        T = zeros(nx,1);
		T(10) = 1;
		b = NaN(nx,1);
		figure(1);clf;
        h_line = plot(T);
        h_title = title('0');
        % ylim([0 1.5]);
        for tidx=1:nt
			for xidx = 1:nx % construct d, the right hand side of the system of linear equations
				delta = s(xidx);
				if xidx ~= 1;  delta = delta + B(M(xidx),M(xidx-1))*(T(xidx-1)-T(xidx)); end
				if xidx ~= nx; delta = delta + B(M(xidx),M(xidx+1))*(T(xidx+1)-T(xidx)); end
				T(xidx) = T(xidx) + delta;
			end
			
			%% Thomas algorithm, sweeps up from 2 to N and then down from N to 1
			% Forward sweep, index 2:
			b(1) = 1 + A1(M(1),M(2));
			w = -A1(M(2),M(1))/b(1);
			b(2) = 1 + A1(M(2),M(3)) + A1(M(2),M(1)) + w*A1(M(1),M(2));
			T(2) = T(2) - w*T(1);
			% Forward sweep, indices 3 to nx-1:
			for xidx = 3:nx-1
				w = -A1(M(xidx),M(xidx-1))/b(xidx-1);
				b(xidx) = 1 + A1(M(xidx),M(xidx-1)) + A1(M(xidx),M(xidx+1));
				b(xidx) = b(xidx) + w*A1(M(xidx-1),M(xidx));
				T(xidx) = T(xidx) - w*T(xidx-1);
			end
			% Forward sweep, index nx:
			w = -A1(M(nx),M(nx-1))/b(nx-1);
			b(nx) = 1 + A1(M(nx),M(nx-1));
			b(nx) = b(nx) + w*A1(M(nx-1),M(nx));
			T(nx) = T(nx) - w*T(nx-1);

			% Back sweep, index nx:
			T(nx) = T(nx)/b(nx);
			% Back sweep, remaining indices:
			for xidx = nx-1:-1:1
				T(xidx) = (T(xidx) + A1(M(xidx),M(xidx+1))*T(xidx+1))/b(xidx);
			end
			h_line.YData = T;
            h_title.String = num2str(tidx);
            drawnow;
        end
        
        
        
        
        
        
        
        
        
        
    case 2
        N = nx*ny;
        M = ones(nx,ny);
%         M(5:15,23:30) = 2;
%         M(1:2,1:2) = 2;
        ABS = zeros(nx,ny); % Power deposited per unit length/area/volume per unit time
%         ABS(16,22) = 100;
        ABS(3,3) = 100;
        
        A1xy_full     = zeros(N); % In xy form
        A2xy_full     = zeros(N); % In xy form
        A2yx_full     = zeros(N); % In yx form
        Bxy_full      = -eye(N); % In xy form
        sxy           = zeros(N,1); % In xy form
        for i=1:N % i is the xy linear index
            [ix,iy] = ind2sub([nx,ny],i);
            j = sub2ind([ny nx],iy,ix); % j is the yx linear index
            heatsinkedvoxel = heatsinkedboundary && (ix == 1 || iy == 1 || ix == nx || iy == ny);
            
            %% x couplings
            if ix ~= 1 && ix ~= nx
                TC_eff_posx = 2*TC(M(i))*TC(M(i+1))/(TC(M(i))+TC(M(i+1)));
                if isnan(TC_eff_posx)
                    TC_eff_posx = 0;
                end
                TC_eff_negx = 2*TC(M(i))*TC(M(i-1))/(TC(M(i))+TC(M(i-1)));
                if isnan(TC_eff_negx)
                    TC_eff_negx = 0;
                end
                A1xy_full(i,i)   = A1xy_full(i,i)   + dt*(TC_eff_posx + TC_eff_negx)/(2*VHC(M(i))*dx^2);
                A1xy_full(i,i-1) = A1xy_full(i,i-1) - dt*TC_eff_negx/(2*VHC(M(i))*dx^2);
                A1xy_full(i,i+1) = A1xy_full(i,i+1) - dt*TC_eff_posx/(2*VHC(M(i))*dx^2);
                Bxy_full(i,i)    = Bxy_full(i,i  )  + dt*(TC_eff_posx + TC_eff_negx)/(2*VHC(M(i))*dx^2);
                Bxy_full(i,i-1)  = Bxy_full(i,i-1)  - dt*TC_eff_negx/(2*VHC(M(i))*dx^2);
                Bxy_full(i,i+1)  = Bxy_full(i,i+1)  - dt*TC_eff_posx/(2*VHC(M(i))*dx^2);
            elseif ~heatsinkedvoxel
                if ix == 1
                    TC_eff_posx = 2*TC(M(i))*TC(M(i+1))/(TC(M(i))+TC(M(i+1)));
                    if isnan(TC_eff_posx)
                        TC_eff_posx = 0;
                    end
                    A1xy_full(i,i)   = A1xy_full(i,i)   + dt*TC_eff_posx/(2*VHC(M(i))*dx^2);
                    A1xy_full(i,i+1) = A1xy_full(i,i+1) - dt*TC_eff_posx/(2*VHC(M(i))*dx^2);
                    Bxy_full(i,i)    = Bxy_full(i,i  )  + dt*TC_eff_posx/(2*VHC(M(i))*dx^2);
                    Bxy_full(i,i+1)  = Bxy_full(i,i+1)  - dt*TC_eff_posx/(2*VHC(M(i))*dx^2);
                else
                    TC_eff_negx = 2*TC(M(i))*TC(M(i-1))/(TC(M(i))+TC(M(i-1)));
                    if isnan(TC_eff_negx)
                        TC_eff_negx = 0;
                    end
                    A1xy_full(i,i)   = A1xy_full(i,i)   + dt*TC_eff_negx/(2*VHC(M(i))*dx^2);
                    A1xy_full(i,i-1) = A1xy_full(i,i-1) - dt*TC_eff_negx/(2*VHC(M(i))*dx^2);
                    Bxy_full(i,i)    = Bxy_full(i,i  )  + dt*TC_eff_negx/(2*VHC(M(i))*dx^2);
                    Bxy_full(i,i-1)  = Bxy_full(i,i-1)  - dt*TC_eff_negx/(2*VHC(M(i))*dx^2);
                end
            end
            
            %% y couplings
            if iy ~= 1 && iy ~= ny
                TC_eff_posy = 2*TC(M(i))*TC(M(i+nx))/(TC(M(i))+TC(M(i+nx)));
                if isnan(TC_eff_posy)
                    TC_eff_posy = 0;
                end
                TC_eff_negy = 2*TC(M(i))*TC(M(i-nx))/(TC(M(i))+TC(M(i-nx)));
                if isnan(TC_eff_negy)
                    TC_eff_negy = 0;
                end
                A2xy_full(i,i)    = A2xy_full(i,i)    + dt*(TC_eff_posy + TC_eff_negy)/(2*VHC(M(i))*dy^2);
                A2xy_full(i,i-nx) = A2xy_full(i,i-nx) - dt*TC_eff_negy /(2*VHC(M(i))*dy^2);
                A2xy_full(i,i+nx) = A2xy_full(i,i+nx) - dt*TC_eff_posy /(2*VHC(M(i))*dy^2);
                A2yx_full(j,j)    = A2yx_full(j,j)    + dt*(TC_eff_posy + TC_eff_negy)/(2*VHC(M(i))*dy^2);
                A2yx_full(j,j-1)  = A2yx_full(j,j-1)  - dt*TC_eff_negy /(2*VHC(M(i))*dy^2);
                A2yx_full(j,j+1)  = A2yx_full(j,j+1)  - dt*TC_eff_posy /(2*VHC(M(i))*dy^2);
                Bxy_full(i,i)     = Bxy_full(i,i   )  + dt*(TC_eff_posy + TC_eff_negy)/(2*VHC(M(i))*dy^2);
                Bxy_full(i,i-nx)  = Bxy_full(i,i-nx)  - dt*TC_eff_negy /(2*VHC(M(i))*dy^2);
                Bxy_full(i,i+nx)  = Bxy_full(i,i+nx)  - dt*TC_eff_posy /(2*VHC(M(i))*dy^2);
            elseif ~heatsinkedvoxel
                if iy == 1
                    TC_eff_posy = 2*TC(M(i))*TC(M(i+nx))/(TC(M(i))+TC(M(i+nx)));
                    if isnan(TC_eff_posy)
                        TC_eff_posy = 0;
                    end
                    A2xy_full(i,i)    = A2xy_full(i,i)    + dt*TC_eff_posy/(2*VHC(M(i))*dy^2);
                    A2xy_full(i,i+nx) = A2xy_full(i,i+nx) - dt*TC_eff_posy/(2*VHC(M(i))*dy^2);
                    A2yx_full(j,j)    = A2yx_full(j,j)    + dt*TC_eff_posy/(2*VHC(M(i))*dy^2);
                    A2yx_full(j,j+1)  = A2yx_full(j,j+1)  - dt*TC_eff_posy/(2*VHC(M(i))*dy^2);
                    Bxy_full(i,i)     = Bxy_full(i,i   )  + dt*TC_eff_posy/(2*VHC(M(i))*dy^2);
                    Bxy_full(i,i+nx)  = Bxy_full(i,i+nx)  - dt*TC_eff_posy/(2*VHC(M(i))*dy^2);
                else
                    TC_eff_negy = 2*TC(M(i))*TC(M(i-nx))/(TC(M(i))+TC(M(i-nx)));
                    if isnan(TC_eff_negy)
                        TC_eff_negy = 0;
                    end
                    A2xy_full(i,i)    = A2xy_full(i,i)    + dt*TC_eff_negy/(2*VHC(M(i))*dy^2);
                    A2xy_full(i,i-nx) = A2xy_full(i,i-nx) - dt*TC_eff_negy/(2*VHC(M(i))*dy^2);
                    A2yx_full(j,j)    = A2yx_full(j,j)    + dt*TC_eff_negy/(2*VHC(M(i))*dy^2);
                    A2yx_full(j,j-1)  = A2yx_full(j,j-1)  - dt*TC_eff_negy/(2*VHC(M(i))*dy^2);
                    Bxy_full(i,i)     = Bxy_full(i,i   )  + dt*TC_eff_negy/(2*VHC(M(i))*dy^2);
                    Bxy_full(i,i-nx)  = Bxy_full(i,i-nx)  - dt*TC_eff_negy/(2*VHC(M(i))*dy^2);
                end
            end
            
            %% Absorption
            if ~heatsinkedvoxel
                sxy(i) = ABS(i)*dt/VHC(M(i));
            end
        end
        A1xy = sparse(A1xy_full);
        A2xy = sparse(A2xy_full);
        A2yx = sparse(A2yx_full);
        Bxy  = sparse(Bxy_full);
        IplusA1xy = eye(N) + A1xy;
        minusBxyminusA2xy = - Bxy - A2xy;
        IplusA2yx = eye(N) + A2yx;
        
        figure(2);clf;
        imagesc(M.');
        axis equal
        axis tight
        axis xy
        
        Txy = zeros(nx*ny,1); % In xy form
        Txy([1 end]) = 0;
        figure(1);clf;
        h_im = imagesc(x,y,reshape(Txy,[nx ny]).');
        colorbar;
        axis equal
        axis tight
        axis xy
        h_title = title('0');
%         caxis([0 0.2]);
        for i=1:nt
            v1xy = IplusA1xy\(minusBxyminusA2xy*Txy + sxy); % This is the intermediate temperature estimate in xy form
            h_im.CData = reshape(v1xy,[nx ny]).';
%             Txy = v1xy;
            rhs2xy = v1xy + A2xy*Txy;
            rhs2yx = reshape(permute(reshape(rhs2xy,[nx ny]),[2 1]),[N 1]);
            Tyx = IplusA2yx\rhs2yx; % New temperature
            h_im.CData = reshape(Tyx,[ny nx]);
            Txy = reshape(permute(reshape(Tyx,[ny nx]),[2 1]),[N 1]);
            colorbar;
            h_title.String = num2str(i);
            drawnow;
        end
        
end
