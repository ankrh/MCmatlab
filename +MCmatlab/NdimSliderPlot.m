function [h_f,h_a] = NdimSliderPlot(data,varargin)
%   arguments
%     data                                    % An N-dimensional array of the values of the dependent variable.
%   end
%   arguments (Repeating)
%     varargin                                % See parameters below.
%   end

  %% Parse optional Name,Value pairs
  p = inputParser;
  addOptional(p,'nFig',NaN,@(x)isScalarPositiveIntegerOrNaN(x)); % An integer denoting which figure number to make the plot in. If NaN, the figure number will be generated automatically to avoid overwriting existing figures.
  addOptional(p,'axisValues',{},@(x)iscell(x)); % A (1,N) cell array, each element containing a 1D array of the values of that independent variable.
  addOptional(p,'axisLabels',{},@(x)iscell(x)); % A (1,N+1) cell array, each element containing the label of one dimension. The first N elements are for the N independent variables and the last element is for the dependent variable.
  addOptional(p,'plotLimits',[NaN NaN],@(x)(numel(x) == 2 && isnumeric(x) && all(~isinf(x))) && (any(isnan(x)) || x(2) > x(1))); % [lower, upper] limits for the data variable. If either limit is NaN, the plotting function will scale automatically to the minimum or maximum value in the data set. Minimum plot limits for log plots are always set to 40 dB under the max.
  addOptional(p,'axisEqual',false,@(x)islogical(x)); % If true, will set the scale of the axis dimensions to be equal. This is useful, for example, when the axis dimensions are spatial (x,y) or (x,y,z) and you want to preserve the actual spatial aspect ratios in the visualization
  addOptional(p,'reversedAxes',[],@isUniqueNaturalNumberArray); % An array of unique natural numbers denoting which of the plot axes to be reversed. For example, [1 3] will reverse the first and third axes, e.g., x and z.
  addOptional(p,'linColormap',MCmatlab.inferno,@(x)(isnumeric(x) && ismatrix(x))); % The colormap to use for linear plot mode
  addOptional(p,'logColormap',MCmatlab.makec2f,@(x)(isnumeric(x) && ismatrix(x))); % The colormap to use for log plot mode
  addOptional(p,'axisDims',[],@isUniqueNaturalNumberArray); % An array of unique natural numbers denoting which dimensions of the data array should be used for axes in the visualization. All other dimensions that are not singleton dimensions are adjusted with sliders.
  addOptional(p,'slicePositions',[],@(x)(isnumeric(x) && all(x>=0) && all(x<=1))); % For geometry illustration plots
  addOptional(p,'docked',true,@(x)(islogical(x))); % Should created figure be docked?
  % The following parameter is used only for plot types used in MCmatlab
  % and is not intended to be for general-purpose use:
  addOptional(p,'indexLabels',{},@(x)iscell(x)); % For geometry illustration plots

  parse(p,varargin{:});

  if isempty(p.Results.axisValues)
    warning('No axis values provided for generalizedPlot.');
    nD = find(size(data) > 1,1,'last');
    for iD = 1:nD
      axisValues{iD} = 1:size(data,iD);
    end
  else
    axisValues = p.Results.axisValues;
  end

  if isempty(p.Results.axisLabels)
    warning('No axis labels provided for generalizedPlot.');
    nD = find(size(data) > 1,1,'last');
    for iD = 1:nD
      axisLabels{iD} = ['Axis ' num2str(iD)];
    end
    axisLabels{iD+1} = 'Data';
  else
    axisLabels = p.Results.axisLabels;
  end

  if isempty(p.Results.slicePositions)
    slicePositions = [0.5 1 0];
    if ismember(2,p.Results.reversedAxes)
      slicePositions(2) = 0;
    end
    if ismember(3,p.Results.reversedAxes)
      slicePositions(3) = 1;
    end
  else
    slicePositions = p.Results.slicePositions;
  end

  %% Convert data to single precision or uint8 to reduce memory footprint
  if isempty(p.Results.indexLabels)
    data = single(data);
  else
    data = uint8(data);
  end

  %% Find out which dimensions are non-singleton, which correspond to axes in the figure, and which are slider dimensions
  if isrow(data) && numel(axisLabels) == 2 && numel(axisValues) == 1
    data = data.'; %% Assume data should be treated as a column vector
  end
  nD = find(size(data) > 1,1,'last'); % Number of dimensions in data (not counting trailing singleton dimensions, because MATLAB automatically cuts them off)
  NonSiDims = find(size(data) > 1); % Indices of non-singleton dimensions
  nNonSi = numel(NonSiDims); % Number of non-singleton dimensions in data
  AxDims = p.Results.axisDims;
  if isempty(AxDims) % If the caller did not specify any axis dimensions, use (up to) the first three non-singleton dimensions
    AxDims = NonSiDims;
    AxDims(4:end) = [];
  end
  nAx = numel(AxDims); % Number of (independent variable) axes to use in the figure

  if nAx > 3
    error('Error: You cannot specify more than three axis dimensions');
  elseif any(setdiff(AxDims,NonSiDims))
    error('Error: You have specified a singleton dimension to be an axis dimension');
  elseif isempty(data)
    error('Error: Input array is empty')
  elseif isscalar(data)
    error('Error: Input array is a scalar')
  elseif numel(axisLabels) < nD + 1
    error('Error: Not enough axis labels specified.');
  elseif numel(axisValues) < nD
    error('Error: Not enough axis value arrays specified');
  end

  if nAx == 3 % If 3D plot
    nSl = nNonSi; % Number of sliders needed
  else % 1D or 2D plot
    nSl = nNonSi - nAx; % Number of sliders needed
  end
  nNonSl = nNonSi - nSl; % Number of non-slider non-singleton dimensions

  %% Permute the data array
  % We will place the axis dimensions first, then the non-singleton
  % dimensions and finally the singleton dimensions, which will
  % subsequently automatically be discarded by MATLAB. The permuted array
  % therefore has no singleton dimensions.
  permuteDimOrder = [AxDims setdiff(NonSiDims,AxDims) setdiff(1:nD,NonSiDims)];
  if ~isscalar(permuteDimOrder)
    data = permute(data,permuteDimOrder);
    axisValues = axisValues(permuteDimOrder(1:nNonSi));
    axisLabels = axisLabels([permuteDimOrder(1:nNonSi) end]);
  end
  
  %% Create figure and put in the necessary uicontrols such as sliders
  if isnan(p.Results.nFig)
    h_f = figure;
  else
    h_f = figure(p.Results.nFig);
  end
  if p.Results.docked
    set(h_f,'WindowStyle','Docked');
  else
    set(h_f,'WindowStyle','normal');
  end
  clf reset;
  h_f.Color = 'w';
  h_a = axes;
  set(h_a,'fontsize',18);

  h_sliders = gobjects(1,nSl);
  h_centerbuttons = gobjects(1,nSl);
  h_slidertexts = gobjects(1,nSl);
  subscriptIdxs = cell(1,nNonSi);
  for iD = 1:nNonSl % For all non-slider dimensions
    subscriptIdxs{iD} = 1:size(data,iD);
  end
  for iSl = 1:nSl % For all slider dimensions
    iD = iSl + nNonSl;
    n = size(data,iD); % Size of the dimension
    if nAx == 3 && iSl <= 3
      n = n + 1; % For 3D visualizations we work with a data array that is padded in the first three axes
    end
    i = round(n/2);
    h_centerbuttons(iSl) = uicontrol('Parent',h_f,'Style','pushbutton','String','Center','Position',[10,iSl*20-10,40,20]);
    warning('off','MATLAB:hg:UIControlSliderStepValueDifference');
    h_sliders(iSl) = uicontrol('Parent',h_f,'Style','slider','Position',[55,iSl*20-10,100,20],...
      'value',i, 'min',1, 'max',n,'SliderStep',[1/(n-1) max(1/(n-1),0.1)]);
    warning('on','MATLAB:hg:UIControlSliderStepValueDifference');
    string = [axisLabels{iD} ' = ' num2str(axisValues{iD}(i),'%.3g')];
    h_slidertexts(iSl) = uicontrol('style','text','String',string,...
      'horizontalAlignment','left','BackgroundColor','w','Position',[160,7 + (iSl-1)*20,5.5*numel(string),20]);
    subscriptIdxs{iD} = i;
  end
  h_checkbox = uicontrol('Parent',h_f,'Style','checkbox','String','log plot','BackgroundColor','w','Position',[10,nSl*20 + 13,55,20]);

  %% Do the plotting
  if ~isnan(p.Results.plotLimits(1)) 
    plotLimLower = p.Results.plotLimits(1);
  else
    plotLimLower = min(data(:));
  end
  if ~isnan(p.Results.plotLimits(2)) 
    plotLimUpper = p.Results.plotLimits(2);
  else
    plotLimUpper = max(data(:));
  end
  if plotLimUpper <= plotLimLower
    plotLimUpper = plotLimLower + 1;
  end

  switch nAx
    case 1 % 1D plot (one independent variable axis and one dependent variable axis)
      h_1D = plot(axisValues{1},data(subscriptIdxs{:}),'linewidth',2);
      grid on; grid minor;
      xlabel(axisLabels{1});
      ylabel(axisLabels{end});

      if p.Results.axisEqual
        h_a.DataAspectRatio = [1 1 1];
      end
      if any(p.Results.reversedAxes == 1)
        h_a.XDir = 'reverse';
      end
      if any(p.Results.reversedAxes == 2)
        h_a.YDir = 'reverse';
      end
      if isempty(p.Results.indexLabels)
        h_a.YLim = [plotLimLower plotLimUpper];
      else
        set(h_checkbox,'Visible','off');
      end

      h_f.UserData = data;

      vars = struct('h_1D',h_1D,...
        'h_centerbuttons',h_centerbuttons,...
        'h_sliders',h_sliders,'h_slidertexts',h_slidertexts, ...
        'h_checkbox',h_checkbox,'linColormap',p.Results.linColormap,...
        'logColormap',p.Results.logColormap,...
        'axisValues',{axisValues},'axisLabels',{axisLabels},...
        'plotLimits',p.Results.plotLimits);
    case 2 % 2D plot (two independent variable axes - the dependent variable is shown with color)
      h_2D = imagesc(axisValues{1},axisValues{2},data(subscriptIdxs{:}).');
      xlabel(axisLabels{1});
      ylabel(axisLabels{2});
      title(axisLabels{end});

      if p.Results.axisEqual
        h_a.DataAspectRatio = [1 1 1];
      end
      if any(p.Results.reversedAxes == 1)
        h_a.XDir = 'reverse';
      end
      if any(p.Results.reversedAxes == 2)
        h_a.YDir = 'reverse';
      else
        h_a.YDir = 'normal';
      end
      if isempty(p.Results.indexLabels)
        colorbar;
        colormap(p.Results.linColormap);
        caxis(h_a,[plotLimLower plotLimUpper]);
      else
        set(h_checkbox,'Visible','off');
        % To get the legend to display properly, we have to cheat MATLAB a
        % bit by introducing some invisible patches with the colors we want:
        nM = numel(p.Results.indexLabels);
        h_patches = gobjects(1,nM);
        colormap(p.Results.linColormap(1:nM,:));
        for iColor = 1:nM
          h_patches(iColor) = patch(NaN,NaN,NaN,iColor);
        end
        legend(h_patches,p.Results.indexLabels,'AutoUpdate','off','FontSize',8,'Location','eastoutside');
      end

      h_f.UserData = data;

      vars = struct('h_2D',h_2D,...
        'h_centerbuttons',h_centerbuttons,...
        'h_sliders',h_sliders,'h_slidertexts',h_slidertexts, ...
        'h_checkbox',h_checkbox,'linColormap',p.Results.linColormap,...
        'logColormap',p.Results.logColormap,...
        'axisValues',{axisValues},'axisLabels',{axisLabels},...
        'plotLimits',p.Results.plotLimits);
    case 3 % 3D plot (three independent variable axes - the dependent variable is shown with color)
      xraw = axisValues{1};
      yraw = axisValues{2};
      zraw = axisValues{3};
      axisValues{1} = [axisValues{1} axisValues{1}(end)];
      axisValues{2} = [axisValues{2} axisValues{2}(end)];
      axisValues{3} = [axisValues{3} axisValues{3}(end)];

      dx = xraw(2)-xraw(1);
      dy = yraw(2)-yraw(1);
      dz = zraw(2)-zraw(1);

      % Voxel corner positions are calculated, since surfaces are drawn based
      % on corner coordinates, not center coordinates. Vectors and matrix are
      % padded accordingly:
      x = round([(xraw - dx/2) , (max(xraw) + dx/2)],15);
      y = round([(yraw - dy/2) , (max(yraw) + dy/2)],15);
      z = round([(zraw - dz/2) , (max(zraw) + dz/2)],15);

      dataSizePadded = size(data);
      dataSizePadded(1:3) = dataSizePadded(1:3) + 1;
      h_f.UserData = zeros(dataSizePadded,'like',data);

      h_f.UserData(1:end-1,1:end-1,1:end-1,:) = data(:,:,:,:);
      h_f.UserData(end,:,:,:) = h_f.UserData(end-1,:,:,:);
      h_f.UserData(:,end,:,:) = h_f.UserData(:,end-1,:,:);
      h_f.UserData(:,:,end,:) = h_f.UserData(:,:,end-1,:);

      xl = x(1); % x low
      xh = x(end); % x high
      yl = y(1); % y low
      yh = y(end); % y high
      zl = z(1); % z low
      zh = z(end); % z high

      [nx,ny,nz,~] = size(h_f.UserData);
      if any(p.Results.reversedAxes == 1)
        h_a.XDir = 'reverse';
        xb = xl;
        xbi = 1;
      else
        xb = xh;
        xbi = nx;
      end
      if any(p.Results.reversedAxes == 2)
        h_a.YDir = 'reverse';
        yb = yl;
        ybi = 1;
      else
        yb = yh;
        ybi = ny;
      end
      if any(p.Results.reversedAxes == 3)
        h_a.ZDir = 'reverse';
        zb = zh;
        zbi = nz;
      else
        zb = zl;
        zbi = 1;
      end

      set(h_sliders(1),'Value',1+(xbi-1)*slicePositions(1));
      set(h_sliders(2),'Value',1+(ybi-1)*slicePositions(2));
      set(h_sliders(3),'Value',1+(zbi-1)*slicePositions(3));
      xsi = floor(get(h_sliders(1),'Value'));
      ysi = floor(get(h_sliders(2),'Value'));
      zsi = floor(get(h_sliders(3),'Value'));
      xs = x(xsi);
      ys = y(ysi);
      zs = z(zsi);

      h_surfxback  = surface(repmat(xb,ny,nz),repmat(y', 1,nz),repmat(z ,ny, 1),squeeze(h_f.UserData(xbi,:,:,subscriptIdxs{4:end})),'LineStyle','none');
      h_surfyback  = surface(repmat(x', 1,nz),repmat(yb,nx,nz),repmat(z ,nx, 1),squeeze(h_f.UserData(:,ybi,:,subscriptIdxs{4:end})),'LineStyle','none');
      h_surfzback  = surface(repmat(x', 1,ny),repmat(y ,nx, 1),repmat(zb,nx,ny),squeeze(h_f.UserData(:,:,zbi,subscriptIdxs{4:end})),'LineStyle','none');
      h_surfxslice = surface(repmat(xs,ny,nz),repmat(y', 1,nz),repmat(z ,ny, 1),squeeze(h_f.UserData(xsi,:,:,subscriptIdxs{4:end})),'LineStyle','none');
      h_surfyslice = surface(repmat(x', 1,nz),repmat(ys,nx,nz),repmat(z ,nx, 1),squeeze(h_f.UserData(:,ysi,:,subscriptIdxs{4:end})),'LineStyle','none');
      h_surfzslice = surface(repmat(x', 1,ny),repmat(y ,nx, 1),repmat(zs,nx,ny),squeeze(h_f.UserData(:,:,zsi,subscriptIdxs{4:end})),'LineStyle','none');

      if isempty(p.Results.indexLabels)
        colorbar;
        linecolor = [0.5 0.5 0.5]; % Assume gray lines around all slices
        caxis(h_a,[plotLimLower plotLimUpper]);
        colormap(p.Results.linColormap);
      else
        set(h_checkbox,'Visible','off');
        linecolor = [0 0 0]; % Black lines around slices
        % To get the legend to display properly, we have to cheat MATLAB a
        % bit by introducing some invisible patches with the colors we want:
        nM = numel(p.Results.indexLabels);
        h_patches = gobjects(1,nM);
        colormap(p.Results.linColormap(1:nM,:));
        for iColor = 1:nM
          h_patches(iColor) = patch(NaN,NaN,NaN,iColor);
        end
        legend(h_patches,p.Results.indexLabels,'AutoUpdate','off','FontSize',8,'Location','eastoutside');
      end

      line([xb xb xb xb xb],[yl yh yh yl yl],[zh zh zl zl zh],'Color',linecolor);
      line([xl xh xh xl xl],[yb yb yb yb yb],[zh zh zl zl zh],'Color',linecolor);
      line([xl xl xh xh xl],[yl yh yh yl yl],[zb zb zb zb zb],'Color',linecolor);
      h_xline = line([xs xs xs xs xs xs xs],[ys yh yh yl yl ys ys],[zl zl zh zh zl zl zh],'Color',linecolor);
      h_yline = line([xl xl xh xh xl xl xh],[ys ys ys ys ys ys ys],[zs zh zh zl zl zs zs],'Color',linecolor);
      h_zline = line([xs xh xh xl xl xs xs],[yl yl yh yh yl yl yh],[zs zs zs zs zs zs zs],'Color',linecolor);

      axis tight
      view(3);
      rotate3d on
      setAxes3DPanAndZoomStyle(zoom(gca),gca,'camera');

      xlabel(axisLabels{1});
      ylabel(axisLabels{2});
      zlabel(axisLabels{3});
      title(axisLabels{end});
      if p.Results.axisEqual
        h_a.DataAspectRatio = [1 1 1];
      end
      drawnow; % To set the automatic CLims
      set(h_a,'CLimMode','manual');

      vars = struct('h_surfxback',h_surfxback,'h_surfyback',h_surfyback,'h_surfzback',h_surfzback,...
        'h_surfxslice',h_surfxslice,'h_surfyslice',h_surfyslice,'h_surfzslice',h_surfzslice,...
        'h_xline',h_xline,'h_yline',h_yline,'h_zline',h_zline,...
        'x',x,'y',y,'z',z,'xbi',xbi,'ybi',ybi,'zbi',zbi,...
        'h_centerbuttons',h_centerbuttons,...
        'h_sliders',h_sliders,'h_slidertexts',h_slidertexts, ...
        'h_checkbox',h_checkbox,'linColormap',p.Results.linColormap,...
        'logColormap',p.Results.logColormap,...
        'axisValues',{axisValues},'axisLabels',{axisLabels},...
        'plotLimits',p.Results.plotLimits);
  end

  set([h_centerbuttons h_sliders h_checkbox],'Callback',{@MCmatlab.NdimSliderPlotRedraw,vars});
  
  for iSl = 1:nSl % For all slider dimensions
    iD = iSl + nNonSl;
    string = [axisLabels{iD} ' = ' num2str(axisValues{iD}(floor(get(h_sliders(iSl),'Value'))),'%.3g')];
    set(h_slidertexts(iSl),'String',string);
    prevpos = get(h_slidertexts(iSl),'Position');
    newpos = [prevpos(1:2) 5.5*numel(string) prevpos(4)];
    set(h_slidertexts(iSl),'Position',newpos);
  end
end

function y = isUniqueNaturalNumberArray(x)
  y = isa(x,'double') && (isempty(x) || isrow(x)) && all(rem(x,1) == 0) && all(x > 0) && numel(x) == numel(unique(x));
end

function y = isScalarPositiveIntegerOrNaN(x)
  y = isscalar(x) && (isnan(x) || (isnumeric(x) && isfinite(x) && floor(x) == x && x > 0));
end
