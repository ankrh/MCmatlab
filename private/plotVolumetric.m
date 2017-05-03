function plotVolumetric(xraw,yraw,zraw,Mraw,varargin)

clf;
h_f = gcf;
if ~isprop(h_f,'M')
    addprop(h_f,'M');
end

if ~isprop(h_f,'Mlog')
addprop(h_f,'Mlog');
end

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

h_f.M = padarray(double(Mraw),[1 1 1],'replicate','post');
h_f.Mlog = log10(h_f.M);

xl = x(1); % x low
xh = x(end); % x high
yl = y(1); % y low
yh = y(end); % y high
zl = z(1); % z low
zh = z(end); % z high

[X,Y,Z] = ndgrid(x,y,z);

warning('off','MATLAB:hg:UIControlSliderStepValueDifference');
h_slider1 = uicontrol('Parent',h_f,'Style','slider','Position',[30,20,200,20],...
              'value',length(x), 'min',1, 'max',length(x),'SliderStep',[1/(length(x)-1) 0.1]);
h_slider2 = uicontrol('Parent',h_f,'Style','slider','Position',[30,40,200,20],...
              'value',length(y), 'min',1, 'max',length(y),'SliderStep',[1/(length(y)-1) 0.1]);
h_slider3 = uicontrol('Parent',h_f,'Style','slider','Position',[30,60,200,20],...
              'value',length(z), 'min',1, 'max',length(z),'SliderStep',[1/(length(z)-1) 0.1]);
warning('on','MATLAB:hg:UIControlSliderStepValueDifference');
uicontrol('style','text','String','x','Position',[10,18,20,20])
uicontrol('style','text','String','y','Position',[10,38,20,20])
uicontrol('style','text','String','z','Position',[10,58,20,20])
h_checkbox1 = uicontrol('Parent',h_f,'Style','checkbox','Position',[70,90,20,20]);
h_checkbox1text = uicontrol('style','text','String','log10 plot','Position',[16,87,50,20]);

h_surfxmax   = surface(squeeze(X(end,:,:)),squeeze(Y(end,:,:)),squeeze(Z(end,:,:)),squeeze(h_f.M(end,:,:)),'LineStyle','none');
h_surfymax   = surface(squeeze(X(:,end,:)),squeeze(Y(:,end,:)),squeeze(Z(:,end,:)),squeeze(h_f.M(:,end,:)),'LineStyle','none');
h_surfzmax   = surface(squeeze(X(:,:,end)),squeeze(Y(:,:,end)),squeeze(Z(:,:,end)),squeeze(h_f.M(:,:,end)),'LineStyle','none');
h_surfxslice = surface(squeeze(X(end,:,:)),squeeze(Y(end,:,:)),squeeze(Z(end,:,:)),squeeze(h_f.M(end,:,:)),'LineStyle','none');
h_surfyslice = surface(squeeze(X(:,end,:)),squeeze(Y(:,end,:)),squeeze(Z(:,end,:)),squeeze(h_f.M(:,end,:)),'LineStyle','none');
h_surfzslice = surface(squeeze(X(:,:,end)),squeeze(Y(:,:,end)),squeeze(Z(:,:,end)),squeeze(h_f.M(:,:,end)),'LineStyle','none');

line([xh xh xh xh xh],[yl yh yh yl yl],[zh zh zl zl zh],'Color','k');
line([xl xh xh xl xl],[yh yh yh yh yh],[zh zh zl zl zh],'Color','k');
line([xl xl xh xh xl],[yl yh yh yl yl],[zh zh zh zh zh],'Color','k');
h_xline = line([xh xh xh xh xh xh xh],[yh yh yh yl yl yh yh],[zl zl zh zh zl zl zh],'Color','k');
h_yline = line([xl xl xh xh xl xl xh],[yh yh yh yh yh yh yh],[zh zh zh zl zl zh zh],'Color','k');
h_zline = line([xh xh xh xl xl xh xh],[yl yl yh yh yl yl yh],[zh zh zh zh zh zh zh],'Color','k');

if ~isempty(varargin)
    tissueList = varargin{1};
    set(h_checkbox1,'Visible','off');
    set(h_checkbox1text,'Visible','off');
    set([h_surfxmax h_surfymax h_surfzmax h_surfxslice h_surfyslice h_surfzslice],'CDataMapping','direct');
    colormap(hsv(length(tissueList)));
    colorbar('TickLabels',{tissueList.name},'Ticks',(1:length(tissueList))+0.5);
else
    colormap(parula(1024));
    colorbar
    if min(h_f.M(:)) ~= max(h_f.M(:))
        caxis([min(h_f.M(:)) max(h_f.M(:))]);
    else
        caxis([min(h_f.M(:)) min(h_f.M(:))+1]);
    end
end
axis tight
axis equal
xlabel('x [cm]')
ylabel('y [cm]')
zlabel('z [cm]')
set(gca,'fontsize',18)
set(gca,'ZDir','reverse')
view(3)
setAxes3DPanAndZoomStyle(zoom(gca),gca,'camera')

vars = struct('h_checkbox1',h_checkbox1,...
    'h_slider1',h_slider1,'h_slider2',h_slider2,'h_slider3',h_slider3,...
    'x',x,'y',y,'z',z,...
    'h_surfxmax',h_surfxmax,'h_surfymax',h_surfymax,'h_surfzmax',h_surfzmax,...
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
            set(vars.h_surfxmax  ,'CData',squeeze(h_f.Mlog(end,:,:)));
            set(vars.h_surfymax  ,'CData',squeeze(h_f.Mlog(:,end,:)));
            set(vars.h_surfzmax  ,'CData',squeeze(h_f.Mlog(:,:,end)));
            set(vars.h_surfxslice,'CData',squeeze(h_f.Mlog(xsi,:,:)));
            set(vars.h_surfyslice,'CData',squeeze(h_f.Mlog(:,ysi,:)));
            set(vars.h_surfzslice,'CData',squeeze(h_f.Mlog(:,:,zsi)));
            if ~isempty(event) % event is empty if callback was made because M was changed. In that case we don't want to renormalize the color scale.
                caxis([max(h_f.Mlog(:))-4 max(h_f.Mlog(:))]);
            end
        else
            set(vars.h_surfxmax  ,'CData',squeeze(h_f.M(end,:,:)));
            set(vars.h_surfymax  ,'CData',squeeze(h_f.M(:,end,:)));
            set(vars.h_surfzmax  ,'CData',squeeze(h_f.M(:,:,end)));
            set(vars.h_surfxslice,'CData',squeeze(h_f.M(xsi,:,:)));
            set(vars.h_surfyslice,'CData',squeeze(h_f.M(:,ysi,:)));
            set(vars.h_surfzslice,'CData',squeeze(h_f.M(:,:,zsi)));
            if ~isempty(event) % event is empty if callback was made because M was changed. In that case we don't want to renormalize the color scale.
                caxis([min(h_f.M(:)) max(h_f.M(:))]);
            end
        end
    case vars.h_slider1
        set(vars.h_surfxslice,'XData',xs*ones(length(vars.y),length(vars.z)));
        if plotLog
            set(vars.h_surfxslice,'CData',squeeze(h_f.Mlog(xsi,:,:)));
        else
            set(vars.h_surfxslice,'CData',squeeze(h_f.M(xsi,:,:)));
        end
        set(vars.h_xline,'XData',xs*ones(1,7));
        set(vars.h_zline,'XData',[xs xh xh xl xl xs xs]);
    case vars.h_slider2
        set(vars.h_surfyslice,'YData',ys*ones(length(vars.x),length(vars.z)));
        if plotLog
            set(vars.h_surfyslice,'CData',squeeze(h_f.Mlog(:,ysi,:)));
        else
            set(vars.h_surfyslice,'CData',squeeze(h_f.M(:,ysi,:)));
        end
        set(vars.h_yline,'YData',ys*ones(1,7));
        set(vars.h_xline,'YData',[ys yh yh yl yl ys ys]);
    case vars.h_slider3
        set(vars.h_surfzslice,'ZData',zs*ones(length(vars.x),length(vars.y)));
        if plotLog
            set(vars.h_surfzslice,'CData',squeeze(h_f.Mlog(:,:,zsi)));
        else
            set(vars.h_surfzslice,'CData',squeeze(h_f.M(:,:,zsi)));
        end
        set(vars.h_zline,'ZData',zs*ones(1,7));
        set(vars.h_yline,'ZData',[zs zh zh zl zl zs zs]);
end

drawnow;

return
end