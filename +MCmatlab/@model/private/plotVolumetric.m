function h_f = plotVolumetric(nFig,xraw,yraw,zraw,Mraw,varargin)

h_f = figure(nFig);
h_f.Color = 'w';

clf;

axes

dx = xraw(2)-xraw(1);
dy = yraw(2)-yraw(1);
dz = zraw(2)-zraw(1);
% Voxel corner positions are calculated, since surfaces are drawn based
% on corner coordinates, not center coordinates. Vectors and matrix are
% padded accordingly:
x = round([(xraw - dx/2) , (max(xraw) + dx/2)],15);
y = round([(yraw - dy/2) , (max(yraw) + dy/2)],15);
z = round([(zraw - dz/2) , (max(zraw) + dz/2)],15);

h_f.UserData = zeros(size(Mraw)+[1 1 1]);
h_f.UserData(1:end-1,1:end-1,1:end-1) = single(Mraw);
h_f.UserData(end,:,:) = h_f.UserData(end-1,:,:);
h_f.UserData(:,end,:) = h_f.UserData(:,end-1,:);
h_f.UserData(:,:,end) = h_f.UserData(:,:,end-1);

% padarray(single(Mraw),[1 1 1],'replicate','post'); % The single data type is used to conserve memory
clear Mraw;
[nx,ny,nz] = size(h_f.UserData);

xl = x(1); % x low
xh = x(end); % x high
yl = y(1); % y low
yh = y(end); % y high
zl = z(1); % z low
zh = z(end); % z high

xyzaxes = false; % Assume axes are not necessarily spatial axes and should not be labeled or set "axis equal"
checkboxvisible = true; % Assume the log10 checkbox should be visible
directmapping = false; % Assume CDataMapping should not be set to direct
reverseZ = false; % Assume z axis is not inverted
fromZero = false; % Assume that the minimum of the color scale should not necessarily be zero
colormap(inferno); % Assume we want to use the inferno colormap
linecolor = [0.5 0.5 0.5]; % Assume gray lines around all slices
if(any(strcmp(varargin,'MCmatlab_GeometryIllustration')))
  mediaProperties = varargin{2};
  nM = length(mediaProperties);
  h_legendobjects = zeros(nM,1);
  for i=1:nM
    h_legendobjects(i) = surface(zeros(2),zeros(2),zeros(2),i*ones(2)); % A hack to get MATLAB to display the legend the way I want
  end
  legend(h_legendobjects,{mediaProperties.name},'FontSize',12,'AutoUpdate','off');%,'Location','northeastoutside');

  xyzaxes = true;
  reverseZ = true;
  checkboxvisible = false;
  directmapping = true;
  colormap(lines(nM));
  linecolor = [0 0 0]; % Black lines around slices
elseif(any(strcmp(varargin,'MCmatlab_fromZero')))
  xyzaxes = true;
  reverseZ = true;
  fromZero = true;
  colorbar;
elseif(any(strcmp(varargin,'MCmatlab_heat')))
  xyzaxes = true;
  reverseZ = true;
  colormap(hot);
  colorbar;
else
  colorbar;
end

slicePositionVarargin = find(strcmpi(varargin,'slicePositions'),1); % Empty array if no slicePositions
if(~isempty(slicePositionVarargin))
  slicePositions = varargin{slicePositionVarargin+1}; % A 1-by-3 array of relative (x,y,z) slice positions (0: slice made at lowest value, 1: slice made at highest value)
  xsi = max(1,min(nx,round((nx-1)*slicePositions(1) + 1)));
  ysi = max(1,min(ny,round((ny-1)*slicePositions(2) + 1)));
  zsi = max(1,min(nz,round((nz-1)*slicePositions(3) + 1)));
elseif reverseZ
  xsi = round(nx/2);
  ysi = ny;
  zsi = nz;
else
  xsi = round(nx/2);
  ysi = ny;
  zsi = 1;
end

xs = x(xsi);
ys = y(ysi);
zs = z(zsi);

if reverseZ
  zb = zh; % "z back"
  zbi = nz; % z back index
else
  zb = zl;
  zbi = 1;
end

if ~directmapping
  if min(h_f.UserData(:)) ~= max(h_f.UserData(:))
    if fromZero
      caxis([0 max(h_f.UserData(:))]);
    else
      caxis([min(h_f.UserData(:)) max(h_f.UserData(:))]);
    end
  else
    caxis([min(h_f.UserData(:)) min(h_f.UserData(:))+1]);
  end
end

warning('off','MATLAB:hg:UIControlSliderStepValueDifference');
h_slider1 = uicontrol('Parent',h_f,'Style','slider','Position',[30,20,200,20],...
                      'value',xsi, 'min',1, 'max',nx,'SliderStep',[1/(nx-1) 0.1]);
h_slider2 = uicontrol('Parent',h_f,'Style','slider','Position',[30,40,200,20],...
                      'value',ysi, 'min',1, 'max',ny,'SliderStep',[1/(ny-1) 0.1]);
h_slider3 = uicontrol('Parent',h_f,'Style','slider','Position',[30,60,200,20],...
                      'value',zsi, 'min',1, 'max',nz,'SliderStep',[1/(nz-1) 0.1]);
warning('on','MATLAB:hg:UIControlSliderStepValueDifference');
uicontrol('style','text','String','x','BackgroundColor','w','Position',[10,18,20,20])
uicontrol('style','text','String','y','BackgroundColor','w','Position',[10,38,20,20])
uicontrol('style','text','String','z','BackgroundColor','w','Position',[10,58,20,20])
h_checkbox1 = uicontrol('Parent',h_f,'Style','checkbox','BackgroundColor','w','Position',[70,90,20,20]);
h_checkbox1text = uicontrol('style','text','String','log10 plot','BackgroundColor','w','Position',[16,87,50,20]);
if ~checkboxvisible
  set(h_checkbox1,'Visible','off');
  set(h_checkbox1text,'Visible','off');
end

h_surfxback  = surface(repmat(xh,ny,nz),repmat(y', 1,nz),repmat(z ,ny, 1),squeeze(h_f.UserData(end,:,:)),'LineStyle','none');
h_surfyback  = surface(repmat(x', 1,nz),repmat(yh,nx,nz),repmat(z ,nx, 1),squeeze(h_f.UserData(:,end,:)),'LineStyle','none');
h_surfzback  = surface(repmat(x', 1,ny),repmat(y ,nx, 1),repmat(zb,nx,ny),squeeze(h_f.UserData(:,:,zbi)),'LineStyle','none');
h_surfxslice = surface(repmat(xs,ny,nz),repmat(y', 1,nz),repmat(z ,ny, 1),squeeze(h_f.UserData(xsi,:,:)),'LineStyle','none');
h_surfyslice = surface(repmat(x', 1,nz),repmat(ys,nx,nz),repmat(z ,nx, 1),squeeze(h_f.UserData(:,ysi,:)),'LineStyle','none');
h_surfzslice = surface(repmat(x', 1,ny),repmat(y ,nx, 1),repmat(zs,nx,ny),squeeze(h_f.UserData(:,:,zsi)),'LineStyle','none');

if directmapping
  set([h_surfxback h_surfyback h_surfzback h_surfxslice h_surfyslice h_surfzslice],'CDataMapping','direct');
  set(gca,'CLimMode','manual');
end

line([xh xh xh xh xh],[yl yh yh yl yl],[zh zh zl zl zh],'Color',linecolor);
line([xl xh xh xl xl],[yh yh yh yh yh],[zh zh zl zl zh],'Color',linecolor);
line([xl xl xh xh xl],[yl yh yh yl yl],[zb zb zb zb zb],'Color',linecolor);
h_xline = line([xs xs xs xs xs xs xs],[ys yh yh yl yl ys ys],[zl zl zh zh zl zl zh],'Color',linecolor);
h_yline = line([xl xl xh xh xl xl xh],[ys ys ys ys ys ys ys],[zs zh zh zl zl zs zs],'Color',linecolor);
h_zline = line([xs xh xh xl xl xs xs],[yl yl yh yh yl yl yh],[zs zs zs zs zs zs zs],'Color',linecolor);

axis tight
if xyzaxes
  axis equal
  xlabel('x [cm]')
  ylabel('y [cm]')
  zlabel('z [cm]')
end
set(gca,'fontsize',18)
if reverseZ
  set(gca,'ZDir','reverse')
end
view(3)

if ~verLessThan('matlab','9.0')
  setAxes3DPanAndZoomStyle(zoom(gca),gca,'camera');
end

vars = struct('h_checkbox1',h_checkbox1,...
  'h_slider1',h_slider1,'h_slider2',h_slider2,'h_slider3',h_slider3,...
  'x',x,'y',y,'z',z,...
  'h_surfxback',h_surfxback,'h_surfyback',h_surfyback,'h_surfzback',h_surfzback,...
  'h_surfxslice',h_surfxslice,'h_surfyslice',h_surfyslice,'h_surfzslice',h_surfzslice,...
  'h_xline',h_xline,'h_yline',h_yline,'h_zline',h_zline,'fromZero',fromZero,'zbi',zbi);

set(h_checkbox1,'Callback',{@redrawVolumetric,vars});
set(h_slider1,'Callback',{@redrawVolumetric,vars});
set(h_slider2,'Callback',{@redrawVolumetric,vars});
set(h_slider3,'Callback',{@redrawVolumetric,vars});

rotate3d on
end

