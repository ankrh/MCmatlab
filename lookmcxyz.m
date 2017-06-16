function lookmcxyz(name)
%% Load data from makeTissue.m and MonteCarlo.m
load(['./Data/' name '.mat']);
load(['./Data/' name '_MCoutput.mat'],'F');

%% Make tissue plot
figure(1);
plotVolumetric(x,y,z,T,tissueList);
title('Tissue type illustration');
drawnow;

%% Make tissue properties plot
figure(7);clf;
plotTissueProperties(tissueList);
drawnow;

%% Make fluence rate plot
figure(2);
plotVolumetric(x,y,z,F);
title('Fluence rate (Intensity) [W/cm^2/W.incident] ')
drawnow;

%% Make power absorption plot
mua_vec = [tissueList.mua];
figure(3);
plotVolumetric(x,y,z,mua_vec(T).*F);
title('Normalized absorbed power per unit volume [W/cm^3/W.incident] ')
drawnow;
return