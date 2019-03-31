function mindtmax = calcdtmax(T,TC,VHC,dx,dy,dz)

TC  = single(TC(T));
effectiveTCx = zeros(size(T) + [1 0 0],'single');
effectiveTCy = zeros(size(T) + [0 1 0],'single');
effectiveTCz = zeros(size(T) + [0 0 1],'single');
effectiveTCx(2:end-1,:,:) = 2*TC(1:end-1,:,:).*TC(2:end,:,:)./(TC(1:end-1,:,:)+TC(2:end,:,:));
effectiveTCy(:,2:end-1,:) = 2*TC(:,1:end-1,:).*TC(:,2:end,:)./(TC(:,1:end-1,:)+TC(:,2:end,:));
effectiveTCz(:,:,2:end-1) = 2*TC(:,:,1:end-1).*TC(:,:,2:end)./(TC(:,:,1:end-1)+TC(:,:,2:end));
clear TC

effectiveTCx(isnan(effectiveTCx)) = 0; % Neighboring insulating voxels would return NaN but should just be 0
effectiveTCy(isnan(effectiveTCy)) = 0;
effectiveTCz(isnan(effectiveTCz)) = 0;

individual_dtmax = VHC(T)./(effectiveTCx(1:end-1,:,:)/dx^2 + effectiveTCx(2:end,:,:)/dx^2 ...
                          + effectiveTCy(:,1:end-1,:)/dy^2 + effectiveTCy(:,2:end,:)/dy^2 ...
                          + effectiveTCz(:,:,1:end-1)/dz^2 + effectiveTCz(:,:,2:end)/dz^2);
mindtmax = double(min(individual_dtmax(:)));
end