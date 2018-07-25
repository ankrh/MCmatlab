function [li , w] = convertTempSensorPos(tempSensorPositions,x,y,z)
% To interpolate the temperature sensor temperatures during heat
% simulation, the C script is going to need temperature values from the 8
% surrounding voxels, calculated in the MEX file with this formula:
% T_sensor = 
%     (1-wx)*(1-wy)*(1-wz)*T(ix  ,iy  ,iz  ) + (1-wx)*(1-wy)*wz*T(ix  ,iy  ,iz+1) +
%     (1-wx)*   wy *(1-wz)*T(ix  ,iy+1,iz  ) + (1-wx)*   wy *wz*T(ix  ,iy+1,iz+1) + 
%        wx *(1-wy)*(1-wz)*T(ix+1,iy  ,iz  ) +    wx *(1-wy)*wz*T(ix+1,iy  ,iz+1) +
%        wx *   wy *(1-wz)*T(ix+1,iy+1,iz  ) +    wx *   wy *wz*T(ix+1,iy+1,iz+1);
% 
% This function calculates the linear index corresponding to the subscript
% indices (ix,iy,iz) and the weights wx, wy and wz. The index is calculated
% in C-notation, that is, starting from 0.

nx = length(x);
ny = length(y);
nz = length(z);
dx = x(2)-x(1);
dy = y(2)-y(1);
dz = z(2)-z(1);

nS = size(tempSensorPositions,1); % Number of sensors
li = zeros(nS,1); % Linear index of interpolation group corner voxel
w = zeros(nS,3); % Weights in the x, y and z directions used in the interpolation

for i=1:nS
    relposx = tempSensorPositions(i,1) - x(1);
    if relposx < 0
        ix = 0;
        w(i,1) = 0;
    elseif relposx < x(end) - x(1)
        ix = floor(relposx/dx);
        w(i,1) = mod(relposx/dx,1);
    else
        ix = nx-2;
        w(i,1) = 1;
    end

    relposy = tempSensorPositions(i,2) - y(1);
    if relposy < 0
        iy = 0;
        w(i,2) = 0;
    elseif relposy < y(end) - y(1)
        iy = floor(relposy/dy);
        w(i,2) = mod(relposy/dy,1);
    else
        iy = ny-2;
        w(i,2) = 1;
    end

    relposz = tempSensorPositions(i,3) - z(1);
    if relposz < 0
        iz = 0;
        w(i,3) = 0;
    elseif relposz < z(end) - z(1)
        iz = floor(relposz/dz);
        w(i,3) = mod(relposz/dz,1);
    else
        iz = nz-2;
        w(i,3) = 1;
    end

    li(i) = ix + ...
            iy*nx + ...
            iz*nx*ny;
end

end