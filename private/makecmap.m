function cmap = makecmap()

f1  = 0.7; % color adjustments
f2  = 0.5; % color adjustments
f3  = 0.3; % color adjustments

cmap(1,:) = [0 0 0]; 
cmap(2,:) = [1 1 1];
cmap(3,:) = [1 0 0];
cmap(4,:) = [1 0.8 0.8];
cmap(5,:) = [0.5 0.2 0.2];
cmap(6,:) = [f1 1 f1]; % skull
cmap(7,:) = [f2 f2 f2]; % gray matter
cmap(8,:) = [0.5 1 1]; % white matter
cmap(9,:) = [1 .5 .5]; % standard tissue

end
