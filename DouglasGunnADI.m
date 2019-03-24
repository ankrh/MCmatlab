calctype = 2; % 1: 1D, 2:2D, 3:3D

insulatingboundary = false;

dt = 0.0001;

dx = 0.1;
dy = 0.05;
dz = 0.08;

TC = [10 ; 1];

VHC = [1 ; 2];

nx = 30;
ny = 40;
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

        A1_full = eye(N);
        B_full = eye(N);
        B_rhs = zeros(N,1);
        for i=1:N
            ix = i;
            insulatedboundaryvoxel = insulatingboundary && (ix == 1 || ix == nx);
            if ix ~= 1 && ix ~= nx
                TC_eff_posx = 2*TC(M(i))*TC(M(i+1))/(TC(M(i))+TC(M(i+1)));
                if isnan(TC_eff_posx)
                    TC_eff_posx = 0;
                end
                TC_eff_negx = 2*TC(M(i))*TC(M(i-1))/(TC(M(i))+TC(M(i-1)));
                if isnan(TC_eff_negx)
                    TC_eff_negx = 0;
                end
                A1_full(i,i)   = A1_full(i,i)   + dt*(TC_eff_posx + TC_eff_negx)/(2*VHC(M(i))*dx^2);
                A1_full(i,i-1) = A1_full(i,i-1) - dt*TC_eff_negx /(2*VHC(M(i))*dx^2);
                A1_full(i,i+1) = A1_full(i,i+1) - dt*TC_eff_posx /(2*VHC(M(i))*dx^2);
                B_full(i,i)   = B_full(i,i)   - dt*(TC_eff_posx + TC_eff_negx)/(2*VHC(M(i))*dx^2);
                B_full(i,i-1) = B_full(i,i-1) + dt*TC_eff_negx /(2*VHC(M(i))*dx^2);
                B_full(i,i+1) = B_full(i,i+1) + dt*TC_eff_posx /(2*VHC(M(i))*dx^2);
            elseif insulatedboundaryvoxel
                if ix == 1
                    TC_eff_posx = 2*TC(M(i))*TC(M(i+1))/(TC(M(i))+TC(M(i+1)));
                    if isnan(TC_eff_posx)
                        TC_eff_posx = 0;
                    end
                    A1_full(i,i)   = A1_full(i,i)   + dt*TC_eff_posx/(2*VHC(M(i))*dx^2);
                    A1_full(i,i+1) = A1_full(i,i+1) - dt*TC_eff_posx /(2*VHC(M(i))*dx^2);
                    B_full(i,i)   = B_full(i,i)   - dt*TC_eff_posx/(2*VHC(M(i))*dx^2);
                    B_full(i,i+1) = B_full(i,i+1) + dt*TC_eff_posx /(2*VHC(M(i))*dx^2);
                else
                    TC_eff_negx = 2*TC(M(i))*TC(M(i-1))/(TC(M(i))+TC(M(i-1)));
                    if isnan(TC_eff_negx)
                        TC_eff_negx = 0;
                    end
                    A1_full(i,i)   = A1_full(i,i)   + dt*TC_eff_negx /(2*VHC(M(i))*dx^2);
                    A1_full(i,i-1) = A1_full(i,i-1) - dt*TC_eff_negx /(2*VHC(M(i))*dx^2);
                    B_full(i,i)   = B_full(i,i)   - dt*TC_eff_negx /(2*VHC(M(i))*dx^2);
                    B_full(i,i-1) = B_full(i,i-1) + dt*TC_eff_negx /(2*VHC(M(i))*dx^2);
                end
            end
            if ~insulatedboundaryvoxel
                B_rhs(i) = ABS(i)*dt/VHC(M(i));
            end
        end
        A1 = sparse(A1_full);
        B = sparse(B_full);

        T = zeros(nx,1);
        T([1 end]) = 0;
        figure(1);clf;
        h_line = plot(T);
        % ylim([0 1.5]);
        for i=1:nt
            T = A1\(B*T + B_rhs);
            h_line.YData = T;
            pause(0.1);
        end
        
        
    case 2
        N = nx*ny;
        M = ones(nx,ny);
        M(5:15,23:30) = 2;
        ABS = zeros(nx,ny); % Power deposited per unit length/area/volume per unit time
        ABS(8,22) = 100;
        
        A1_full     = zeros(N);
        A2_full     = zeros(N);
        B_full      = eye(N);
        s           = zeros(N,1);
        for i=1:N
            [ix,iy] = ind2sub([nx,ny],i);
            insulatedboundaryvoxel = insulatingboundary && (ix == 1 || iy == 1 || ix == nx || iy == ny);
            
            %% x couplings
            if ix ~= 1 && ix ~= nx
                TC_eff_posx = 2*TC(M(i))*TC(M(i+1))/(TC(M(i))+TC(M(i+1)));
                TC_eff_negx = 2*TC(M(i))*TC(M(i-1))/(TC(M(i))+TC(M(i-1)));
                A1_full(i,i)   = A1_full(i,i)   + dt*(TC_eff_posx + TC_eff_negx)/(2*VHC(M(i))*dx^2);
                A1_full(i,i-1) = A1_full(i,i-1) - dt*TC_eff_negx /(2*VHC(M(i))*dx^2);
                A1_full(i,i+1) = A1_full(i,i+1) - dt*TC_eff_posx /(2*VHC(M(i))*dx^2);
                B_full(i,i)    = B_full(i,i  )  - dt*(TC_eff_posx + TC_eff_negx)/(2*VHC(M(i))*dx^2);
                B_full(i,i-1)  = B_full(i,i-1)  + dt*TC_eff_negx /(2*VHC(M(i))*dx^2);
                B_full(i,i+1)  = B_full(i,i+1)  + dt*TC_eff_posx /(2*VHC(M(i))*dx^2);
            elseif insulatedboundaryvoxel
                if ix == 1
                    TC_eff_posx = 2*TC(M(i))*TC(M(i+1))/(TC(M(i))+TC(M(i+1)));
                    A1_full(i,i)   = A1_full(i,i)   + dt*TC_eff_posx/(2*VHC(M(i))*dx^2);
                    A1_full(i,i+1) = A1_full(i,i+1) - dt*TC_eff_posx /(2*VHC(M(i))*dx^2);
                    B_full(i,i)    = B_full(i,i  )  - dt*(TC_eff_posx + TC_eff_negx)/(2*VHC(M(i))*dx^2);
                    B_full(i,i+1)  = B_full(i,i+1)  + dt*TC_eff_posx /(2*VHC(M(i))*dx^2);
                else
                    TC_eff_negx = 2*TC(M(i))*TC(M(i-1))/(TC(M(i))+TC(M(i-1)));
                    A1_full(i,i)   = A1_full(i,i)   + dt*TC_eff_negx /(2*VHC(M(i))*dx^2);
                    A1_full(i,i-1) = A1_full(i,i-1) - dt*TC_eff_negx /(2*VHC(M(i))*dx^2);
                    B_full(i,i)    = B_full(i,i  )  - dt*(TC_eff_posx + TC_eff_negx)/(2*VHC(M(i))*dx^2);
                    B_full(i,i-1)  = B_full(i,i-1)  + dt*TC_eff_negx /(2*VHC(M(i))*dx^2);
                end
            end
            
            %% y couplings
            if iy ~= 1 && iy ~= ny
                TC_eff_posy = 2*TC(M(i))*TC(M(i+nx))/(TC(M(i))+TC(M(i+nx)));
                TC_eff_negy = 2*TC(M(i))*TC(M(i-nx))/(TC(M(i))+TC(M(i-nx)));
                A2_full(i,i)    = A2_full(i,i)    + dt*(TC_eff_posy + TC_eff_negy)/(2*VHC(M(i))*dy^2);
                A2_full(i,i-nx) = A2_full(i,i-nx) - dt*TC_eff_negy /(2*VHC(M(i))*dy^2);
                A2_full(i,i+nx) = A2_full(i,i+nx) - dt*TC_eff_posy /(2*VHC(M(i))*dy^2);
                B_full(i,i)     = B_full(i,i   )  - dt*(TC_eff_posy + TC_eff_negy)/(2*VHC(M(i))*dy^2);
                B_full(i,i-nx)  = B_full(i,i-nx)  + dt*TC_eff_negy /(2*VHC(M(i))*dy^2);
                B_full(i,i+nx)  = B_full(i,i+nx)  + dt*TC_eff_posy /(2*VHC(M(i))*dy^2);

                Ay_lhs_full(i,i)    = Ay_lhs_full(i,i)    + dt*(TC_eff_posy + TC_eff_negy)/(2*VHC(M(i))*dy^2);
                Ay_lhs_full(i,i-1)  = Ay_lhs_full(i,i-1)  - dt*TC_eff_negy /(2*VHC(M(i))*dy^2); % To be used after permuting x and y
                Ay_lhs_full(i,i+1)  = Ay_lhs_full(i,i+1)  - dt*TC_eff_posy /(2*VHC(M(i))*dy^2); % To be used after permuting x and y
            elseif insulatedboundaryvoxel
                if iy == 1
                    TC_eff_posy = 2*TC(M(i))*TC(M(i+1))/(TC(M(i))+TC(M(i+1)));
                    A2_full(i,i)    = A2_full(i,i)    + dt*TC_eff_posy/(2*VHC(M(i))*dy^2);
                    A2_full(i,i+nx) = A2_full(i,i+nx) - dt*TC_eff_posy/(2*VHC(M(i))*dy^2);
                    B_full(i,i)     = B_full(i,i   )  - dt*TC_eff_posy/(2*VHC(M(i))*dy^2);
                    B_full(i,i+nx)  = B_full(i,i+nx)  + dt*TC_eff_posy/(2*VHC(M(i))*dy^2);

                    Ay_lhs_full(i,i)    = Ay_lhs_full(i,i)    + dt*TC_eff_posy/(2*VHC(M(i))*dy^2);
                    Ay_lhs_full(i,i+nx) = Ay_lhs_full(i,i+nx) - dt*TC_eff_posy /(2*VHC(M(i))*dy^2);
                else
                    TC_eff_negy = 2*TC(M(i))*TC(M(i-1))/(TC(M(i))+TC(M(i-1)));
                    A2_full(i,i)    = A2_full(i,i)    + dt*TC_eff_negy/(2*VHC(M(i))*dy^2);
                    A2_full(i,i-nx) = A2_full(i,i-nx) - dt*TC_eff_negy/(2*VHC(M(i))*dy^2);
                    B_full(i,i)     = B_full(i,i   )  - dt*TC_eff_negy/(2*VHC(M(i))*dy^2);
                    B_full(i,i-nx)  = B_full(i,i-nx)  + dt*TC_eff_negy/(2*VHC(M(i))*dy^2);

                    Ay_lhs_full(i,i)    = Ay_lhs_full(i,i)    + dt*TC_eff_negy/(2*VHC(M(i))*dy^2);
                    Ay_lhs_full(i,i-nx) = Ay_lhs_full(i,i-nx) - dt*TC_eff_negy /(2*VHC(M(i))*dy^2);
                end
            end
            
            %% Absorption
            if ~insulatedboundaryvoxel
                s(i) = ABS(i)*dt/VHC(M(i));
            end
        end
        A1 = sparse(A1_full);
        A2 = sparse(A2_full);
        B  = sparse(B_full);
        step1_lhs = eye(N) + A1;
        step1_rhs = - B - A2;
        step2_lhs = eye(N) + A2;
        step2_rhs = A2;
        
        figure(2);clf;
        imagesc(M.');
        axis equal
        axis tight
        axis xy
        
        T = zeros(nx*ny,1);
        T([1 end]) = 0;
        figure(1);clf;
        h_im = imagesc(x,y,reshape(T,[nx ny]).');
        colorbar;
        axis equal
        axis tight
        axis xy
%         caxis([0 0.2]);
        for i=1:nt
            v1 = A1\(B*T + s); % This is the intermediate temperature estimate
            h_im.CData = reshape(T,[nx ny]).';
            colorbar;
            drawnow;
        end
        
end
