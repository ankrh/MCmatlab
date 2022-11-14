function [h_f,h_a] = plotFarField(nFig,axisValues,data,axisLabels)
%   arguments
%     nFig       (1,1) double {mustBeInteger} % An integer denoting which figure number to make the plot in
%     axisValues (1,:) cell                   % A (1,N) cell array, each element containing a 1D array of the values of the independent variable, excluding the far field angles theta and phi.
%     data                                    % An N+2 dimensional array of the values of the dependent variable, where the first two dimensions correspond to theta and phi angles.
%     axisLabels (1,:) cell                   % A (1,N) cell array, each element containing the label of one dimension.
%   end

  %% Convert data to single precision to reduce memory footprint
  data = single(data);

  %% Find out which dimensions are non-singleton, which correspond to axes in the figure, and which are slider dimensions
  nD = find(size(data) > 1,1,'last'); % Number of dimensions in data (not counting trailing singleton dimensions, because MATLAB automatically cuts them off)
  NonSiDims = find(size(data) > 1); % Indices of non-singleton dimensions
  nNonSi = numel(NonSiDims); % Number of non-singleton dimensions in data
  AxDims = [1 2];
  nAx = numel(AxDims); % Number of (independent variable) axes to use in the figure

  if isempty(data)
    error('Error: Input array is empty')
  elseif isscalar(data)
    error('Error: Input array is a scalar')
  elseif isvector(data)
    error('Error: Input array is 1D');
  elseif numel(axisLabels) < nD - 2
    error('Error: Not enough axis labels specified.');
  elseif numel(axisValues) < nD - 2
    error('Error: Not enough axis value arrays specified');
  end

  nSl = nNonSi - nAx; % Number of sliders needed
  nNonSl = nNonSi - nSl; % Number of non-slider non-singleton dimensions

  %% Permute the data array
  % We will place the axis dimensions first, then the non-singleton
  % dimensions and finally the singleton dimensions, which will
  % subsequently automatically be discarded by MATLAB. The permuted array
  % therefore has no singleton dimensions.
  permuteDimOrder = [AxDims setdiff(NonSiDims,AxDims) setdiff(1:nD,NonSiDims)];
  data = permute(data,permuteDimOrder);
  axisValues = axisValues(permuteDimOrder(3:nNonSi)-2);
  axisLabels = axisLabels(permuteDimOrder(3:nNonSi)-2);

  %% Create figure and put in the necessary uicontrols such as sliders
  h_f = figure(nFig);
  set(h_f,'WindowStyle','Docked');
  clf reset;
  h_f.Color = 'w';
  h_a = axes;

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
    i = round(n/2);
    h_centerbuttons(iSl) = uicontrol('Parent',h_f,'Style','pushbutton','String','Center','Position',[10,iSl*20-10,40,20]);
    warning('off','MATLAB:hg:UIControlSliderStepValueDifference');
    h_sliders(iSl) = uicontrol('Parent',h_f,'Style','slider','Position',[55,iSl*20-10,100,20],...
      'value',i, 'min',1, 'max',n,'SliderStep',[1/(n-1) max(1/(n-1),0.1)]);
    warning('on','MATLAB:hg:UIControlSliderStepValueDifference');
    string = [axisLabels{iSl} ' = ' num2str(axisValues{iSl}(i),'%.3g')];
    h_slidertexts(iSl) = uicontrol('style','text','String',string,...
      'horizontalAlignment','left','BackgroundColor','w','Position',[160,7 + (iSl-1)*20,5.5*numel(string),20]);
    subscriptIdxs{iD} = i;
  end
  h_checkbox = uicontrol('Parent',h_f,'Style','checkbox','String','log plot','BackgroundColor','w','Position',[10,nSl*20 + 13,55,20]);

  %% Do the plotting
  nFF = size(data,1);
  theta_vec = linspace(0,pi,nFF+1).'; % MCorFMC.farFieldTheta contains the theta values at the centers of the far field pixels, but we will not use those here since we need the corner positions to calculate the solid angle of the pixels
  solidAngle_vec = 2*pi*(cos(theta_vec(1:end-1)) - cos(theta_vec(2:end)))/nFF; % Solid angle extended by each far field pixel, as function of theta
  [ux,uy,minus_uz] = sphere(nFF); % This MATLAB function gives us the correct ux,uy,uz coordinates of the corners of the far field pixels, except that uz has the wrong sign
  uz = -minus_uz;
  sizes = num2cell(size(data));
  data = data./repmat(solidAngle_vec,1,sizes{2:end});
  data(1  ,:,:) = repelem(mean(data(1  ,:,:)),nFF,1);
  data(end,:,:) = repelem(mean(data(end,:,:)),nFF,1);
  maxelement = max(data(:));
  if maxelement == 0
    maxelement = 1;
  end
  h_f.UserData = data;
  
  h_2D = surf(ux,uy,uz,h_f.UserData(subscriptIdxs{:}),'EdgeColor','none');
  xlabel('u_x');
  ylabel('u_y');
  zlabel('u_z');
  h_a.DataAspectRatio = [1 1 1];
  h_a.ZDir = 'reverse';
  colorbar;
  caxis([0 maxelement]);

  linColormap = inferno;
  logColormap = makec2f;
  colormap(linColormap);
  rotate3d on
  setAxes3DPanAndZoomStyle(zoom(gca),gca,'camera');
  set(gca,'fontsize',18);

  vars = struct('h_2D',h_2D,...
    'h_centerbuttons',h_centerbuttons,...
    'h_sliders',h_sliders,'h_slidertexts',h_slidertexts, ...
    'h_checkbox',h_checkbox,'linColormap',linColormap,...
    'logColormap',logColormap,'axisValues',{axisValues},...
    'axisLabels',{axisLabels},'maxelement',maxelement);

  set([h_centerbuttons h_sliders h_checkbox],'Callback',{@redrawFarField,vars});
end