function individual_dtmax = calcdtmax(Tissue,TC,HC,dx,dy,dz)

TC  = TC(Tissue);
effectiveTCx = 2*TC(1:end-1,:,:).*TC(2:end,:,:)./(TC(1:end-1,:,:)+TC(2:end,:,:));
effectiveTCy = 2*TC(:,1:end-1,:).*TC(:,2:end,:)./(TC(:,1:end-1,:)+TC(:,2:end,:));
effectiveTCz = 2*TC(:,:,1:end-1).*TC(:,:,2:end)./(TC(:,:,1:end-1)+TC(:,:,2:end));
clear TC

effectiveTCx(isnan(effectiveTCx)) = 0; % Neighboring insulating voxels would return NaN but should just be 0
effectiveTCy(isnan(effectiveTCy)) = 0;
effectiveTCz(isnan(effectiveTCz)) = 0;

HC = HC(Tissue);
individual_dtmax = HC./(padarray(effectiveTCx,[1 0 0],0,'pre')/dx*dy*dz + padarray(effectiveTCx,[1 0 0],0,'post')/dx*dy*dz ...
                      + padarray(effectiveTCy,[0 1 0],0,'pre')*dx/dy*dz + padarray(effectiveTCy,[0 1 0],0,'post')*dx/dy*dz ...
                      + padarray(effectiveTCz,[0 0 1],0,'pre')*dx*dy/dz + padarray(effectiveTCz,[0 0 1],0,'post')*dx*dy/dz);
end