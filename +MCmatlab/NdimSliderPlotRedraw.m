function generalizedPlotRedraw(src,~,v)
  h_f = src.Parent;
  nSl = numel(v.h_sliders);
  nNonSi = find(size(h_f.UserData) > 1,1,'last');
  nNonSl = nNonSi - nSl;
  if nNonSl == 0
    nAx = 3;
  else
    nAx = nNonSl;
  end

  %% Set slider if source was one of the center buttons
  for iSl = find(src == v.h_centerbuttons)
    iD = iSl + nNonSl;
    set(v.h_sliders(iSl),'Value',(numel(v.axisValues{iD})-1)/2 + 1);
  end
  
  %% Set text if source was one of the sliders or center buttons
  for iSl = [find(src == v.h_sliders) find(src == v.h_centerbuttons)]
    iD = iSl + nNonSl;
    string = [v.axisLabels{iD} ' = ' num2str(v.axisValues{iD}(floor(get(v.h_sliders(iSl),'Value'))),'%.3g')];
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
  switch nAx
    case 1
      h_ax = v.h_1D.Parent;
      v.h_1D.YData = h_f.UserData(subscriptIdxs{:});
      if v.h_checkbox.Value
        h_ax.YScale = 'log';
        ylim(h_ax,[v.maxelement/10000 v.maxelement]);
      else
        h_ax.YScale = 'lin';
        ylim(h_ax,[v.minelement v.maxelement]);
      end
    case 2
      h_ax = v.h_2D.Parent;
      if isgraphics(v.h_2D,'Image')
        v.h_2D.CData = h_f.UserData(subscriptIdxs{:}).';
      else
        v.h_2D.CData = h_f.UserData(subscriptIdxs{:});
      end
      if v.h_checkbox.Visible
        if v.h_checkbox.Value
          h_ax.ColorScale = 'log';
          caxis(h_ax,[v.maxelement/10000 v.maxelement]);
          colormap(v.logColormap);
        else
          h_ax.ColorScale = 'lin';
          caxis(h_ax,[v.minelement v.maxelement]);
          colormap(v.linColormap);
        end
      end
    case 3
      h_ax = v.h_surfxslice.Parent;
      [nx,ny,nz,~] = size(h_f.UserData);
      xl = v.x(1); % x low
      xh = v.x(end); % x high
      yl = v.y(1); % y low
      yh = v.y(end); % y high
      zl = v.z(1); % z low
      zh = v.z(end); % z high
      [xsi,ysi,zsi] = subscriptIdxs{1:3};
      xs = v.x(xsi); ys = v.y(ysi); zs = v.z(zsi);
      if ismember(src,[v.h_sliders(1) v.h_centerbuttons(1)])
        set(v.h_surfxslice,'XData',xs*ones(ny,nz));
        set(v.h_xline,'XData',[xs xs xs xs xs xs xs]);
        set(v.h_zline,'XData',[xs xh xh xl xl xs xs]);
      end
      if ismember(src,[v.h_sliders(2) v.h_centerbuttons(2)])
        set(v.h_surfyslice,'YData',ys*ones(nx,nz));
        set(v.h_yline,'YData',[ys ys ys ys ys ys ys]);
        set(v.h_xline,'YData',[ys yh yh yl yl ys ys]);
      end
      if ismember(src,[v.h_sliders(3) v.h_centerbuttons(3)])
        set(v.h_surfzslice,'ZData',zs*ones(nx,ny));
        set(v.h_zline,'ZData',[zs zs zs zs zs zs zs]);
        set(v.h_yline,'ZData',[zs zh zh zl zl zs zs]);
      end
      if ismember(src,[v.h_sliders([1 4:end]) v.h_centerbuttons([1 4:end])])
        set(v.h_surfxslice,'CData',squeeze(h_f.UserData(xsi,:,:,subscriptIdxs{4:end})));
      end
      if ismember(src,[v.h_sliders([2 4:end]) v.h_centerbuttons([2 4:end])])
        set(v.h_surfyslice,'CData',squeeze(h_f.UserData(:,ysi,:,subscriptIdxs{4:end})));
      end
      if ismember(src,[v.h_sliders([3 4:end]) v.h_centerbuttons([3 4:end])])
        set(v.h_surfzslice,'CData',squeeze(h_f.UserData(:,:,zsi,subscriptIdxs{4:end})));
      end
      if ismember(src,[v.h_sliders(4:end) v.h_centerbuttons(4:end)])
        set(v.h_surfxback,'CData',squeeze(h_f.UserData(v.xbi,:,:,subscriptIdxs{4:end})));
        set(v.h_surfyback,'CData',squeeze(h_f.UserData(:,v.ybi,:,subscriptIdxs{4:end})));
        set(v.h_surfzback,'CData',squeeze(h_f.UserData(:,:,v.zbi,subscriptIdxs{4:end})));
      end
      if v.h_checkbox.Visible
        if v.h_checkbox.Value
          h_ax.ColorScale = 'log';
          caxis(h_ax,[v.maxelement/10000 v.maxelement]);
          colormap(v.logColormap);
        else
          h_ax.ColorScale = 'lin';
          caxis(h_ax,[v.minelement v.maxelement]);
          colormap(v.linColormap);
        end
      end
  end
  drawnow;
end
