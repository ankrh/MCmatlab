function plotVolumetric(x,y,z,M,varargin)
h_f = gcf;
h_slider1 = uicontrol('Parent',h_f,'Style','slider','Position',[30,20,200,20],...
              'value',max(x), 'min',min(x), 'max',max(x));
h_slider2 = uicontrol('Parent',h_f,'Style','slider','Position',[30,40,200,20],...
              'value',max(y), 'min',min(y), 'max',max(y));
h_slider3 = uicontrol('Parent',h_f,'Style','slider','Position',[30,60,200,20],...
              'value',max(z), 'min',min(z), 'max',max(z));
uicontrol('style','text','String','x','Position',[10,18,20,20])
uicontrol('style','text','String','y','Position',[10,38,20,20])
uicontrol('style','text','String','z','Position',[10,58,20,20])
if isempty(varargin)
    h_checkbox1 = uicontrol('Parent',h_f,'Style','checkbox','Position',[70,90,20,20]);
    uicontrol('style','text','String','log10 plot','Position',[16,87,50,20])
    vars = struct('h_checkbox1',h_checkbox1,'h_slider1',h_slider1,'h_slider2',h_slider2,'h_slider3',h_slider3,'x',x,'y',y,'z',z,'M',M);
    set(h_checkbox1,'Callback',{@slider_callback,vars});
else
    vars = struct('h_slider1',h_slider1,'h_slider2',h_slider2,'h_slider3',h_slider3,'x',x,'y',y,'z',z,'M',M,'tissueList',varargin);
end
set(h_slider1,'Callback',{@slider_callback,vars});
set(h_slider2,'Callback',{@slider_callback,vars});
set(h_slider3,'Callback',{@slider_callback,vars});
slider_callback(0,0,vars)
return

function slider_callback(~,~,vars)
if isfield(vars,'h_checkbox1')
    plotLog = get(vars.h_checkbox1,'Value');
    saturationcutoff = prctile(double(vars.M(:)),99.5);
else
    plotLog = 0;
end
xslice = get(vars.h_slider1,'Value');
yslice = get(vars.h_slider2,'Value');
zslice = get(vars.h_slider3,'Value');
x = vars.x;
y = vars.y;
z = vars.z;
h_a = gca;
previoustitle = h_a.Title.String;

if isfield(vars,'tissueList')
    slice(x,y,z,double(vars.M),[xslice max(x)],[yslice max(y)],[zslice max(z)],'nearest');
    tissueList = vars.tissueList;
    % cmap = makecmap();
    % colormap(cmap)
    colormap(hsv(9))
    set(get(gca,'Children'),'LineStyle','none','CDataMapping','direct')
    colorbar('TickLabels',{tissueList.name},'Ticks',(1:length(tissueList))+0.5);
elseif plotLog
    logcutoff = -1;
    slice(x,y,z,max(logcutoff,log10(min(saturationcutoff,double(vars.M)))),[xslice max(x)],[yslice max(y)],[zslice max(z)],'nearest');
    set(get(gca,'Children'),'LineStyle','none','CDataMapping','scaled')
    colorbar
else
    slice(x,y,z,min(saturationcutoff,double(vars.M)),[xslice max(x)],[yslice max(y)],[zslice max(z)],'nearest');
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

line([xslice xslice xslice xslice xslice],[min(y) max(y) max(y) min(y) min(y)],[max(z) max(z) min(z) min(z) max(z)],'Color','k')
line([max(x) max(x) max(x) max(x) max(x)],[min(y) max(y) max(y) min(y) min(y)],[max(z) max(z) min(z) min(z) max(z)],'Color','k')
line([min(x) max(x) max(x) min(x) min(x)],[yslice yslice yslice yslice yslice],[max(z) max(z) min(z) min(z) max(z)],'Color','k')
line([min(x) max(x) max(x) min(x) min(x)],[max(y) max(y) max(y) max(y) max(y)],[max(z) max(z) min(z) min(z) max(z)],'Color','k')
line([min(x) min(x) max(x) max(x) min(x)],[min(y) max(y) max(y) min(y) min(y)],[zslice zslice zslice zslice zslice],'Color','k')
line([min(x) min(x) max(x) max(x) min(x)],[min(y) max(y) max(y) min(y) min(y)],[max(z) max(z) max(z) max(z) max(z)],'Color','k')

line([min(x) max(x)],[yslice yslice],[zslice zslice],'Color','k')
line([xslice xslice],[min(y) max(y)],[zslice zslice],'Color','k')
line([xslice xslice],[yslice yslice],[min(z) max(z)],'Color','k')

set(gca,'fontsize',18)

return
