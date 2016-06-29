function cmap = makecmap(Nt)
% function cmap = makecmap(Nt)
%   Currently, set for makecmap(Nt), where Nt <= 11.
%   You can add more colors as you add tissues to makeTissueList.m.

cmap = zeros(64,3);

dj = 0.05;
for i=1:64
    j = round((i-dj)/64*(Nt-1));
    if      j<=1-dj, cmap(i,:) = [0 0 0]; % escape
    elseif  j<=2-dj, cmap(i,:) = [1 1 1]; % air
    elseif  j<=3-dj, cmap(i,:) = [1 0 0]; % blood
    elseif  j<=4-dj, cmap(i,:) = [1 0.8 0.8]; % dermis
    elseif  j<=5-dj, cmap(i,:) = [0.5 0.2 0.2]; % epidermis
    elseif  j<=6-dj, cmap(i,:) = [0.7 1 0.7]; % skull
    elseif  j<=7-dj, cmap(i,:) = [0.5 0.5 0.5]; % gray matter
    elseif  j<=8-dj, cmap(i,:) = [0.5 1 1]; % white matter
    elseif  j<=9-dj, cmap(i,:) = [0.3 .2 .2]; % hair
    elseif  j<=10-dj, cmap(i,:) = [0 0 1]; % water
    elseif  j<=11-dj, cmap(i,:) = [1 .5 .5]; % standard tissue
    end
end
