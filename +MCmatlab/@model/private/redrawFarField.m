function redrawFarField(src,~,v)
  h_f = gcbf;
  nSl = numel(v.h_sliders);
  nNonSi = find(size(h_f.UserData) > 1,1,'last');
  nNonSl = nNonSi - nSl;

  %% Set slider if source was one of the center buttons
  for iSl = find(src == v.h_centerbuttons)
    set(v.h_sliders(iSl),'Value',(numel(v.axisValues{iSl})-1)/2 + 1);
  end

  %% Set text if source was one of the sliders or center buttons
  for iSl = [find(src == v.h_sliders) find(src == v.h_centerbuttons)]
    string = [v.axisLabels{iSl} ' = ' num2str(v.axisValues{iSl}(floor(get(v.h_sliders(iSl),'Value'))),'%.3g')];
    set(v.h_slidertexts(iSl),'String',string);
    prevpos = get(v.h_slidertexts(iSl),'Position');
    newpos = [prevpos(1:2) 5.5*numel(string) prevpos(4)];
    set(v.h_slidertexts(iSl),'Position',newpos);
  end

  %% Calculate new set of subscript indices to plot
  subscriptIdxs = cell(1,nNonSi);
  for iD = 1:nNonSl % For all non-slider dimensions
    subscriptIdxs{iD} = 1:size(h_f.UserData,iD);
  end
  for iSl = 1:nSl % For all slider dimensions
    iD = iSl + nNonSl;
    subscriptIdxs{iD} = floor(get(v.h_sliders(iSl),'Value'));
  end

  %% Update plot
  h_ax = v.h_2D.Parent;
  v.h_2D.CData = h_f.UserData(subscriptIdxs{:});
  if v.h_checkbox.Value
    h_ax.ColorScale = 'log';
    caxis(h_ax,[v.maxelement/10000 v.maxelement]);
    colormap(v.logColormap);
  else
    h_ax.ColorScale = 'lin';
    caxis(h_ax,[0 v.maxelement]);
    colormap(v.linColormap);
  end

  drawnow;
end
