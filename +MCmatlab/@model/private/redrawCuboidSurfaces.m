function redrawCuboidSurfaces(src,~,v)
lambdatext = [955 ' [nm]']; % The 955 is the unicode number for the lambda character. Because the editors of MATLAB versions prior to 2020a do not support unicode, we have to make the lambda in this sneaky way.

if src == v.h_centerbutton
  set(v.h_slider,'Value',(numel(v.l)-1)/2 + 1);
end

plotLog = get(v.h_checkbox,'Value');
lsi = floor(get(v.h_slider,'Value')); % lambda slice index
ls = v.l(lsi);
string = [lambdatext ' = ' num2str(ls,'%.3g')];

v.h_slidertext.String = string;
prevpos = get(v.h_slidertext,'Position');
newpos = [prevpos(1:2) 5.5*numel(string) prevpos(4)];
set(v.h_slidertext,'Position',newpos);

h_f = get(v.h_checkbox,'Parent');

if plotLog
  colormap(makec2f);
  maxelement = log10(max(max(max([h_f.UserData{:}]))));
  if(isfinite(maxelement))
    caxisvec = [maxelement-4 maxelement];
  else
    caxisvec = [-3 1];
  end
  if ishandle(v.h_xneg)
    v.h_xneg.CData = squeeze(log10(h_f.UserData{1}(:,:,lsi)));
    v.h_xneg.Parent.CLim = caxisvec;
  end
  if ishandle(v.h_yneg)
    v.h_yneg.CData = squeeze(log10(h_f.UserData{2}(:,:,lsi)));
    v.h_yneg.Parent.CLim = caxisvec;
  end
  if ishandle(v.h_zneg)
    v.h_zneg.CData = squeeze(log10(h_f.UserData{3}(:,:,lsi)));
    v.h_zneg.Parent.CLim = caxisvec;
  end
  if ishandle(v.h_xpos)
    v.h_xpos.CData = squeeze(log10(h_f.UserData{4}(:,:,lsi)));
    v.h_xpos.Parent.CLim = caxisvec;
  end
  if ishandle(v.h_ypos)
    v.h_ypos.CData = squeeze(log10(h_f.UserData{5}(:,:,lsi)));
    v.h_ypos.Parent.CLim = caxisvec;
  end
  if ishandle(v.h_zpos)
    v.h_zpos.CData = squeeze(log10(h_f.UserData{6}(:,:,lsi)));
    v.h_zpos.Parent.CLim = caxisvec;
  end
else
  colormap(v.linearColormap);
  maxelement = max(max(max([h_f.UserData{:}])));
  caxisvec = [0 maxelement];
  if ishandle(v.h_xneg)
    v.h_xneg.CData = squeeze(h_f.UserData{1}(:,:,lsi));
    v.h_xneg.Parent.CLim = caxisvec;
  end
  if ishandle(v.h_yneg)
    v.h_yneg.CData = squeeze(h_f.UserData{2}(:,:,lsi));
    v.h_yneg.Parent.CLim = caxisvec;
  end
  if ishandle(v.h_zneg)
    v.h_zneg.CData = squeeze(h_f.UserData{3}(:,:,lsi));
    v.h_zneg.Parent.CLim = caxisvec;
  end
  if ishandle(v.h_xpos)
    v.h_xpos.CData = squeeze(h_f.UserData{4}(:,:,lsi));
    v.h_xpos.Parent.CLim = caxisvec;
  end
  if ishandle(v.h_ypos)
    v.h_ypos.CData = squeeze(h_f.UserData{5}(:,:,lsi));
    v.h_ypos.Parent.CLim = caxisvec;
  end
  if ishandle(v.h_zpos)
    v.h_zpos.CData = squeeze(h_f.UserData{6}(:,:,lsi));
    v.h_zpos.Parent.CLim = caxisvec;
  end
end

drawnow;
return
