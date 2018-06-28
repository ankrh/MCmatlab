function map = GPBGYRcolormap(n)
% GPBGYR: Grey-Purple-Blue-Green-Yellow-Red

if nargin < 1
   n = 1024;
end

values = [
   -1 64 64 64
   0 40 0 60
   256 125 0 115
   512 220 0 180
   768 80 0 255
   1024 0 35 205
   1280 0 125 235
   1536 0 230 255
   1792 0 240 0
   2048 0 200 0
   2304 90 220 0
   2560 233 250 0
   2816 255 244 0
   3072 255 177 0
   3328 255 128 0
   3584 255 96 0
   3840 255 24 0
   4031 255 0 0
   ];

map = interp1(values(:,1), values(:,2:4), linspace(min(values(:,1)),max(values(:,1)),n), 'linear')/255;
