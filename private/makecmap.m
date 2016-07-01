function cmap = makecmap(Nt)
% function cmap = makecmap(Nt)
%   You can add more colors as you add tissues to makeTissueList.m.

cmap = zeros(Nt,3);

for i=1:Nt
    if      i<=1, cmap(i,:) = [0 0 0]; % escape
    elseif  i<=2, cmap(i,:) = [1 1 1]; % air
    elseif  i<=3, cmap(i,:) = [1 0 0]; % blood
    elseif  i<=4, cmap(i,:) = [1 0.8 0.8]; % dermis
    elseif  i<=5, cmap(i,:) = [0.5 0.2 0.2]; % epidermis
    elseif  i<=6, cmap(i,:) = [0.7 1 0.7]; % skull
    elseif  i<=7, cmap(i,:) = [0.5 0.5 0.5]; % gray matter
    elseif  i<=8, cmap(i,:) = [0.5 1 1]; % white matter
    elseif  i<=9, cmap(i,:) = [0.3 .2 .2]; % hair
    elseif  i<=10, cmap(i,:) = [0 0 1]; % water
    elseif  i<=11, cmap(i,:) = [1 .5 .5]; % standard tissue
    end
end
