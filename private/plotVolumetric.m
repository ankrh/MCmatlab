function plotVolumetric(x,y,z,M,varargin)
dx = x(2)-x(1);
dy = y(2)-y(1);
dz = z(2)-z(1);
xpadded = round([(x - dx/2) , (max(x) + dx/2)],15);
ypadded = round([(y - dy/2) , (max(y) + dy/2)],15);
zpadded = round([(z - dz/2) , (max(z) + dz/2)],15);
Mpadded = M;
Mpadded(end+1,:,:) = Mpadded(end,:,:);
Mpadded(:,end+1,:) = Mpadded(:,end,:);
Mpadded(:,:,end+1) = Mpadded(:,:,end);

h_f = gcf;
h_slider1 = uicontrol('Parent',h_f,'Style','slider','Position',[30,20,200,20],...
              'value',max(xpadded), 'min',min(xpadded), 'max',max(xpadded));
h_slider2 = uicontrol('Parent',h_f,'Style','slider','Position',[30,40,200,20],...
              'value',max(ypadded), 'min',min(ypadded), 'max',max(ypadded));
h_slider3 = uicontrol('Parent',h_f,'Style','slider','Position',[30,60,200,20],...
              'value',max(zpadded), 'min',min(zpadded), 'max',max(zpadded));
uicontrol('style','text','String','x','Position',[10,18,20,20])
uicontrol('style','text','String','y','Position',[10,38,20,20])
uicontrol('style','text','String','z','Position',[10,58,20,20])
if isempty(varargin)
    h_checkbox1 = uicontrol('Parent',h_f,'Style','checkbox','Position',[70,90,20,20]);
    uicontrol('style','text','String','log10 plot','Position',[16,87,50,20])
    vars = struct('h_checkbox1',h_checkbox1,'h_slider1',h_slider1,'h_slider2',h_slider2,'h_slider3',h_slider3,'xpadded',xpadded,'ypadded',ypadded,'zpadded',zpadded,'Mpadded',Mpadded);
    set(h_checkbox1,'Callback',{@slider_callback,vars});
else
    vars = struct('h_slider1',h_slider1,'h_slider2',h_slider2,'h_slider3',h_slider3,'xpadded',xpadded,'ypadded',ypadded,'zpadded',zpadded,'Mpadded',Mpadded,'tissueList',varargin);
end
set(h_slider1,'Callback',{@slider_callback,vars});
set(h_slider2,'Callback',{@slider_callback,vars});
set(h_slider3,'Callback',{@slider_callback,vars});
slider_callback(0,0,vars)
return

function slider_callback(~,~,vars)
if isfield(vars,'h_checkbox1')
    plotLog = get(vars.h_checkbox1,'Value');
    saturationcutoff = prctile(double(vars.Mpadded(:)),99.5);
else
    plotLog = 0;
end
xslice = get(vars.h_slider1,'Value');
yslice = get(vars.h_slider2,'Value');
zslice = get(vars.h_slider3,'Value');
xpadded = vars.xpadded;
ypadded = vars.ypadded;
zpadded = vars.zpadded;

h_a = gca;
previoustitle = h_a.Title.String;

if isfield(vars,'tissueList')
    slice(xpadded,ypadded,zpadded,permute(double(vars.Mpadded),[2 1 3]),[xslice max(xpadded)],[yslice max(ypadded)],[zslice max(zpadded)],'nearest');
    tissueList = vars.tissueList;
    % cmap = makecmap();
    % colormap(cmap)
    colormap(hsv(9))
    set(get(gca,'Children'),'LineStyle','none','CDataMapping','direct')
    colorbar('TickLabels',{tissueList.name},'Ticks',(1:length(tissueList))+0.5);
elseif plotLog
    logcutoff = -1;
    slice(xpadded,ypadded,zpadded,max(logcutoff,log10(min(saturationcutoff,permute(double(vars.Mpadded),[2 1 3])))),[xslice max(xpadded)],[yslice max(ypadded)],[zslice max(zpadded)],'nearest');
    set(get(gca,'Children'),'LineStyle','none','CDataMapping','scaled')
    colorbar
else
    slice(xpadded,ypadded,zpadded,min(saturationcutoff,permute(double(vars.Mpadded),[2 1 3])),[xslice max(xpadded)],[yslice max(ypadded)],[zslice max(zpadded)],'nearest');
    set(get(gca,'Children'),'LineStyle','none','CDataMapping','scaled')
    colorbar
end

title(previoustitle);
axis tight
axis equal
xlabel('x [cm]')
ylabel('y [cm]')
zlabel('z [cm]')
set(gca,'ZDir','reverse')

line([xslice xslice xslice xslice xslice],[min(ypadded) max(ypadded) max(ypadded) min(ypadded) min(ypadded)],[max(zpadded) max(zpadded) min(zpadded) min(zpadded) max(zpadded)],'Color','k')
line([max(xpadded) max(xpadded) max(xpadded) max(xpadded) max(xpadded)],[min(ypadded) max(ypadded) max(ypadded) min(ypadded) min(ypadded)],[max(zpadded) max(zpadded) min(zpadded) min(zpadded) max(zpadded)],'Color','k')
line([min(xpadded) max(xpadded) max(xpadded) min(xpadded) min(xpadded)],[yslice yslice yslice yslice yslice],[max(zpadded) max(zpadded) min(zpadded) min(zpadded) max(zpadded)],'Color','k')
line([min(xpadded) max(xpadded) max(xpadded) min(xpadded) min(xpadded)],[max(ypadded) max(ypadded) max(ypadded) max(ypadded) max(ypadded)],[max(zpadded) max(zpadded) min(zpadded) min(zpadded) max(zpadded)],'Color','k')
line([min(xpadded) min(xpadded) max(xpadded) max(xpadded) min(xpadded)],[min(ypadded) max(ypadded) max(ypadded) min(ypadded) min(ypadded)],[zslice zslice zslice zslice zslice],'Color','k')
line([min(xpadded) min(xpadded) max(xpadded) max(xpadded) min(xpadded)],[min(ypadded) max(ypadded) max(ypadded) min(ypadded) min(ypadded)],[max(zpadded) max(zpadded) max(zpadded) max(zpadded) max(zpadded)],'Color','k')

line([min(xpadded) max(xpadded)],[yslice yslice],[zslice zslice],'Color','k')
line([xslice xslice],[min(ypadded) max(ypadded)],[zslice zslice],'Color','k')
line([xslice xslice],[yslice yslice],[min(zpadded) max(zpadded)],'Color','k')

set(gca,'fontsize',18)

return
