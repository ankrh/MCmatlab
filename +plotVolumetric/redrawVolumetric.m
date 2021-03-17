function redrawVolumetric(src,event,vars)

plotLog = get(vars.h_checkbox1,'Value');
xsi = floor(get(vars.h_slider1,'Value')); % x slice index
xs = vars.x(xsi);
ysi = floor(get(vars.h_slider2,'Value')); % y slice index
ys = vars.y(ysi);
zsi = floor(get(vars.h_slider3,'Value')); % z slice index
zs = vars.z(zsi);

xl = vars.x(1); % x low
xh = vars.x(end); % x high
xbi = vars.xbi; % z back index
yl = vars.y(1); % y low
yh = vars.y(end); % y high
zl = vars.z(1); % z low
zh = vars.z(end); % z high
zbi = vars.zbi; % z back index

h_f = get(vars.h_checkbox1,'Parent');

switch src
  case vars.h_checkbox1
    if plotLog
      set(vars.h_surfxback ,'CData',squeeze(log10(h_f.UserData(xbi,:,:))));
      set(vars.h_surfyback ,'CData',squeeze(log10(h_f.UserData(:,end,:))));
      set(vars.h_surfzback ,'CData',squeeze(log10(h_f.UserData(:,:,zbi))));
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
      set(vars.h_surfxback ,'CData',squeeze(h_f.UserData(xbi,:,:)));
      set(vars.h_surfyback ,'CData',squeeze(h_f.UserData(:,end,:)));
      set(vars.h_surfzback ,'CData',squeeze(h_f.UserData(:,:,zbi)));
      set(vars.h_surfxslice,'CData',squeeze(h_f.UserData(xsi,:,:)));
      set(vars.h_surfyslice,'CData',squeeze(h_f.UserData(:,ysi,:)));
      set(vars.h_surfzslice,'CData',squeeze(h_f.UserData(:,:,zsi)));
      if ~isempty(event) % event is empty if callback was made because UserData was changed. In that case we don't want to renormalize the color scale.
        if vars.fromZero
          caxis([0 max(h_f.UserData(:))]);
        else
          caxis([min(h_f.UserData(:)) max(h_f.UserData(:))]);
        end
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
