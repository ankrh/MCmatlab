function plotVolumetric(xraw,yraw,zraw,Mraw,varargin)

clf;
h_f = gcf;

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

h_f.UserData = padarray(single(Mraw),[1 1 1],'replicate','post'); % The single data type is used to conserve memory
clear Mraw;
[nx,ny,nz] = size(h_f.UserData);

xl = x(1); % x low
xh = x(end); % x high
yl = y(1); % y low
yh = y(end); % y high
zl = z(1); % z low
zh = z(end); % z high

checkboxvisible = true; % Assume the log10 checkbox should be visible
directmapping = false; % Assume CDataMapping should not be set to direct
reverseZ = false; % Assume z axis is not inverted
colormap(GPBGYRcolormap); % Assume we want to use the GPBGYR (Grey-Purple-Blue-Green-Yellow-Red) colormap
colorbar;
if ~isempty(varargin)
    plottype = varargin{1};
    if strcmp(plottype,'reverseZTissueIllustration')
        reverseZ = true;
        tissueList = varargin{2};
        checkboxvisible = false;
        directmapping = true;
        colormap(lines(length(tissueList)));
        colorbar('TickLabels',{tissueList.name},'Ticks',(1:length(tissueList))+0.5);
    elseif strcmp(plottype,'reverseZ')
        reverseZ = true;
    end
end

slider3start = 1; % Assume that slider 3 should start at the minimum value
zback = zl; % Assume the back surface in the z direction is the one at high z
if reverseZ
    slider3start = nz;
    zback = zh;
end

if ~directmapping
    if min(h_f.UserData(:)) ~= max(h_f.UserData(:))
        caxis([min(h_f.UserData(:)) max(h_f.UserData(:))]);
    else
        caxis([min(h_f.UserData(:)) min(h_f.UserData(:))+1]);
    end
end

warning('off','MATLAB:hg:UIControlSliderStepValueDifference');
h_slider1 = uicontrol('Parent',h_f,'Style','slider','Position',[30,20,200,20],...
              'value',nx, 'min',1, 'max',nx,'SliderStep',[1/(nx-1) 0.1]);
h_slider2 = uicontrol('Parent',h_f,'Style','slider','Position',[30,40,200,20],...
              'value',ny, 'min',1, 'max',ny,'SliderStep',[1/(ny-1) 0.1]);
h_slider3 = uicontrol('Parent',h_f,'Style','slider','Position',[30,60,200,20],...
              'value',slider3start, 'min',1, 'max',nz,'SliderStep',[1/(nz-1) 0.1]);
warning('on','MATLAB:hg:UIControlSliderStepValueDifference');
uicontrol('style','text','String','x','Position',[10,18,20,20])
uicontrol('style','text','String','y','Position',[10,38,20,20])
uicontrol('style','text','String','z','Position',[10,58,20,20])
h_checkbox1 = uicontrol('Parent',h_f,'Style','checkbox','Position',[70,90,20,20]);
h_checkbox1text = uicontrol('style','text','String','log10 plot','Position',[16,87,50,20]);
if ~checkboxvisible
    set(h_checkbox1,'Visible','off');
    set(h_checkbox1text,'Visible','off');
end

h_surfxback  = surface(repmat(xh,ny,nz),repmat(y', 1,nz),repmat(z    ,ny, 1),squeeze(h_f.UserData(end,:,:)),'LineStyle','none');
h_surfyback  = surface(repmat(x', 1,nz),repmat(yh,nx,nz),repmat(z    ,nx, 1),squeeze(h_f.UserData(:,end,:)),'LineStyle','none');
h_surfzback  = surface(repmat(x', 1,ny),repmat(y ,nx, 1),repmat(zback,nx,ny),squeeze(h_f.UserData(:,:,end)),'LineStyle','none');
h_surfxslice = surface(repmat(xh,ny,nz),repmat(y', 1,nz),repmat(z    ,ny, 1),squeeze(h_f.UserData(end,:,:)),'LineStyle','none');
h_surfyslice = surface(repmat(x', 1,nz),repmat(yh,nx,nz),repmat(z    ,nx, 1),squeeze(h_f.UserData(:,end,:)),'LineStyle','none');
h_surfzslice = surface(repmat(x', 1,ny),repmat(y ,nx, 1),repmat(zback,nx,ny),squeeze(h_f.UserData(:,:,end)),'LineStyle','none');

if directmapping
    set([h_surfxback h_surfyback h_surfzback h_surfxslice h_surfyslice h_surfzslice],'CDataMapping','direct');
end

line([xh xh xh xh xh],[yl yh yh yl yl],[zh zh zl zl zh],'Color','k');
line([xl xh xh xl xl],[yh yh yh yh yh],[zh zh zl zl zh],'Color','k');
line([xl xl xh xh xl],[yl yh yh yl yl],[zh zh zh zh zh],'Color','k');
h_xline = line([xh xh xh xh xh xh xh],[yh yh yh yl yl yh yh],[zl zl zh zh zl zl zh],'Color','k');
h_yline = line([xl xl xh xh xl xl xh],[yh yh yh yh yh yh yh],[zh zh zh zl zl zh zh],'Color','k');
h_zline = line([xh xh xh xl xl xh xh],[yl yl yh yh yl yl yh],[zh zh zh zh zh zh zh],'Color','k');

axis tight
axis equal
xlabel('x [cm]')
ylabel('y [cm]')
zlabel('z [cm]')
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
    'h_xline',h_xline,'h_yline',h_yline,'h_zline',h_zline);

set(h_checkbox1,'Callback',{@callback,vars});
set(h_slider1,'Callback',{@callback,vars});
set(h_slider2,'Callback',{@callback,vars});
set(h_slider3,'Callback',{@callback,vars});
return
end

function callback(src,event,vars)

plotLog = get(vars.h_checkbox1,'Value');
xsi = floor(get(vars.h_slider1,'Value')); % x slice index
xs = vars.x(xsi);
ysi = floor(get(vars.h_slider2,'Value')); % y slice index
ys = vars.y(ysi);
zsi = floor(get(vars.h_slider3,'Value')); % z slice index
zs = vars.z(zsi);

xl = vars.x(1); % x low
xh = vars.x(end); % x high
yl = vars.y(1); % y low
yh = vars.y(end); % y high
zl = vars.z(1); % z low
zh = vars.z(end); % z high

h_f = get(vars.h_checkbox1,'Parent');

switch src
    case vars.h_checkbox1
        if plotLog
            set(vars.h_surfxback ,'CData',squeeze(log10(h_f.UserData(end,:,:))));
            set(vars.h_surfyback ,'CData',squeeze(log10(h_f.UserData(:,end,:))));
            set(vars.h_surfzback ,'CData',squeeze(log10(h_f.UserData(:,:,end))));
            set(vars.h_surfxslice,'CData',squeeze(log10(h_f.UserData(xsi,:,:))));
            set(vars.h_surfyslice,'CData',squeeze(log10(h_f.UserData(:,ysi,:))));
            set(vars.h_surfzslice,'CData',squeeze(log10(h_f.UserData(:,:,zsi))));
            if ~isempty(event) % event is empty if callback was made because UserData was changed. In that case we don't want to renormalize the color scale.
                maxelement = log10(max(h_f.UserData(:)));
                if(isfinite(maxelement))
                    caxis([maxelement-4 maxelement]);
                else
                    caxis([-3 1]);
                end
            end
        else
            set(vars.h_surfxback ,'CData',squeeze(h_f.UserData(end,:,:)));
            set(vars.h_surfyback ,'CData',squeeze(h_f.UserData(:,end,:)));
            set(vars.h_surfzback ,'CData',squeeze(h_f.UserData(:,:,end)));
            set(vars.h_surfxslice,'CData',squeeze(h_f.UserData(xsi,:,:)));
            set(vars.h_surfyslice,'CData',squeeze(h_f.UserData(:,ysi,:)));
            set(vars.h_surfzslice,'CData',squeeze(h_f.UserData(:,:,zsi)));
            if ~isempty(event) % event is empty if callback was made because UserData was changed. In that case we don't want to renormalize the color scale.
                caxis([min(h_f.UserData(:)) max(h_f.UserData(:))]);
            end
        end
    case vars.h_slider1
        set(vars.h_surfxslice,'XData',xs*ones(length(vars.y),length(vars.z)));
        if plotLog
            set(vars.h_surfxslice,'CData',squeeze(log10(h_f.UserData(xsi,:,:))));
        else
            set(vars.h_surfxslice,'CData',squeeze(h_f.UserData(xsi,:,:)));
        end
        set(vars.h_xline,'XData',xs*ones(1,7));
        set(vars.h_zline,'XData',[xs xh xh xl xl xs xs]);
    case vars.h_slider2
        set(vars.h_surfyslice,'YData',ys*ones(length(vars.x),length(vars.z)));
        if plotLog
            set(vars.h_surfyslice,'CData',squeeze(log10(h_f.UserData(:,ysi,:))));
        else
            set(vars.h_surfyslice,'CData',squeeze(h_f.UserData(:,ysi,:)));
        end
        set(vars.h_yline,'YData',ys*ones(1,7));
        set(vars.h_xline,'YData',[ys yh yh yl yl ys ys]);
    case vars.h_slider3
        set(vars.h_surfzslice,'ZData',zs*ones(length(vars.x),length(vars.y)));
        if plotLog
            set(vars.h_surfzslice,'CData',squeeze(log10(h_f.UserData(:,:,zsi))));
        else
            set(vars.h_surfzslice,'CData',squeeze(h_f.UserData(:,:,zsi)));
        end
        set(vars.h_zline,'ZData',zs*ones(1,7));
        set(vars.h_yline,'ZData',[zs zh zh zl zl zs zs]);
end

drawnow;

return
end