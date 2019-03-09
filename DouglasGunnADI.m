insulatingboundary = false;

dt = 0.01*rand;

dx = rand;
dy = 0.09;
dz = 0.08;

TC = [10 10];

VHC = [1 1];

nx = 30;
nt = 1000;
M   = ones(1,nx).';
ABS = zeros(1,nx).';
ABS(20) = 1;

%% Calculate matrices used for Douglas-Gunn ADI method
TCmat = [0 TC(M) 0].';
VHCmat = VHC(M).';

effectiveTCx_mat = 2*TCmat(1:end-1,:,:).*TCmat(2:end,:,:)./(TCmat(1:end-1,:,:)+TCmat(2:end,:,:));
effectiveTCy_mat = 2*TCmat(:,1:end-1,:).*TCmat(:,2:end,:)./(TCmat(:,1:end-1,:)+TCmat(:,2:end,:));
effectiveTCz_mat = 2*TCmat(:,:,1:end-1).*TCmat(:,:,2:end)./(TCmat(:,:,1:end-1)+TCmat(:,:,2:end));
effectiveTCx_mat(isnan(effectiveTCx_mat)) = 0; % Neighboring insulating voxels would return NaN but should just be 0
effectiveTCy_mat(isnan(effectiveTCy_mat)) = 0;
effectiveTCz_mat(isnan(effectiveTCz_mat)) = 0;

Ax_lhs = sparse(numel(M),numel(M));
Ax_lhs(nx+1:nx+1:end) = -effectiveTCx_mat(2:end-1)/(2*dx^2);
Ax_lhs(   1:nx+1:end) = VHCmat/dt + (effectiveTCx_mat(1:end-1) + effectiveTCx_mat(2:end))/(2*dx^2);
Ax_lhs(   2:nx+1:end) = -effectiveTCx_mat(2:end-1)/(2*dx^2);

Ax_rhs = sparse(numel(M),numel(M));
Ax_rhs(nx+1:nx+1:end) = effectiveTCx_mat(2:end-1)/(2*dx^2);
Ax_rhs(   1:nx+1:end) = VHCmat/dt - (effectiveTCx_mat(1:end-1) + effectiveTCx_mat(2:end))/(2*dx^2);
Ax_rhs(   2:nx+1:end) = effectiveTCx_mat(2:end-1)/(2*dx^2);

B_rhs = ABS;

if ~insulatingboundary
    Ax_lhs([1 nx+1 end-nx end]) = [1 0 0 1];
    Ax_rhs([1 nx+1 end-nx end]) = [1 0 0 1];
    B_rhs([1 end])              = [0 0];
end

T = rand(nx,1);
T([1 end]) = 0;
figure(1);clf;
h_line = plot(T);
for i=1:nt
    T = Ax_lhs\(Ax_rhs*T + B_rhs);
    h_line.YData = T;
    pause(0.01);
end

