calctype = 2; % 1: 1D, 2:2D, 3:3D

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
x = dx*(-(nx-1)/2:(nx-1)/2);
y = dy*(-(ny-1)/2:(ny-1)/2);
switch calctype
    case 1
        N = nx;
        M   = ones(nx,1);
        M(12:13) = 2;
        ABS = zeros(nx,1); % Power deposited per unit length/area/volume per unit time
        ABS(10) = 100;

        A1xy_full = eye(N);
        Bxy_full = eye(N);
        B_rhs = zeros(N,1);
        for i=1:N
            ix = i;
            heatsinkedvoxel = heatsinkedboundary && (ix == 1 || ix == nx);
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
                A1xy_full(i,i-1) = A1xy_full(i,i-1) - dt*TC_eff_negx /(2*VHC(M(i))*dx^2);
                A1xy_full(i,i+1) = A1xy_full(i,i+1) - dt*TC_eff_posx /(2*VHC(M(i))*dx^2);
                Bxy_full(i,i)   = Bxy_full(i,i)   - dt*(TC_eff_posx + TC_eff_negx)/(2*VHC(M(i))*dx^2);
                Bxy_full(i,i-1) = Bxy_full(i,i-1) + dt*TC_eff_negx /(2*VHC(M(i))*dx^2);
                Bxy_full(i,i+1) = Bxy_full(i,i+1) + dt*TC_eff_posx /(2*VHC(M(i))*dx^2);
            elseif ~heatsinkedvoxel
                if ix == 1
                    TC_eff_posx = 2*TC(M(i))*TC(M(i+1))/(TC(M(i))+TC(M(i+1)));
                    if isnan(TC_eff_posx)
                        TC_eff_posx = 0;
                    end
                    A1xy_full(i,i)   = A1xy_full(i,i)   + dt*TC_eff_posx/(2*VHC(M(i))*dx^2);
                    A1xy_full(i,i+1) = A1xy_full(i,i+1) - dt*TC_eff_posx /(2*VHC(M(i))*dx^2);
                    Bxy_full(i,i)   = Bxy_full(i,i)   - dt*TC_eff_posx/(2*VHC(M(i))*dx^2);
                    Bxy_full(i,i+1) = Bxy_full(i,i+1) + dt*TC_eff_posx /(2*VHC(M(i))*dx^2);
                else
                    TC_eff_negx = 2*TC(M(i))*TC(M(i-1))/(TC(M(i))+TC(M(i-1)));
                    if isnan(TC_eff_negx)
                        TC_eff_negx = 0;
                    end
                    A1xy_full(i,i)   = A1xy_full(i,i)   + dt*TC_eff_negx /(2*VHC(M(i))*dx^2);
                    A1xy_full(i,i-1) = A1xy_full(i,i-1) - dt*TC_eff_negx /(2*VHC(M(i))*dx^2);
                    Bxy_full(i,i)   = Bxy_full(i,i)   - dt*TC_eff_negx /(2*VHC(M(i))*dx^2);
                    Bxy_full(i,i-1) = Bxy_full(i,i-1) + dt*TC_eff_negx /(2*VHC(M(i))*dx^2);
                end
            end
            if ~heatsinkedvoxel
                B_rhs(i) = ABS(i)*dt/VHC(M(i));
            end
        end
        A1xy = sparse(A1xy_full);
        Bxy = sparse(Bxy_full);

        Txy = zeros(nx,1);
        Txy([1 end]) = 0;
        figure(1);clf;
        h_line = plot(Txy);
        h_title = title('0');
        % ylim([0 1.5]);
        for i=1:nt
            Txy = A1xy\(Bxy*Txy + B_rhs);
            h_line.YData = Txy;
            h_title.String = num2str(i);
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
