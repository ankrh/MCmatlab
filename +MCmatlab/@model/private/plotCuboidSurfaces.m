function h_f = plotCuboidSurfaces(nFig,model,simFluorescence)
  lambdatext = [955 ' [nm]']; % The 955 is the unicode number for the lambda character. Because the editors of MATLAB versions prior to 2020a do not support unicode, we have to make the lambda in this sneaky way.
  
  x = model.G.x;
  y = model.G.y;
  z = model.G.z;

  h_f = figure(nFig);
  set(h_f,'WindowStyle','Docked');
  clf reset;
  h_f.Color = 'w';
  axes

  lincolormap = inferno;

  dx = x(2)-x(1);
  dy = y(2)-y(1);
  dz = z(2)-z(1);

  if simFluorescence
    MCorFMC = model.FMC;
    P_in = dx*dy*dz*sum(model.FMC.sourceDistribution(:));
    fluorescenceOrIncident = 'fluorescence ';
    fluorescenceOrNothing = 'fluorescence ';
  else
    MCorFMC = model.MC;
    P_in = 1;
    fluorescenceOrIncident = 'incident ';
    fluorescenceOrNothing = '';
  end
  [nx,ny,nl] = size(MCorFMC.NI_zneg); % nx and ny are the numbers of bins in the extended box, if using boundaryType 2
  [~ ,nz,~ ] = size(MCorFMC.NI_xpos);
  lsi = round(nl/2); % Lambda "slice" index
  Lx = numel(x)*dx;
  Ly = numel(y)*dy;
  Lz = nz*dz;
  l = MCorFMC.wavelength;

  x = round([(x - dx/2) , (max(x) + dx/2)],15); % Padding
  y = round([(y - dy/2) , (max(y) + dy/2)],15); % Padding
  z = round([(z - dz/2) , (max(z) + dz/2)],15); % Padding
  xl = x(1); % x low
  xh = x(end); % x high
  yl = y(1); % y low
  yh = y(end); % y high
  zl = z(1); % z low
  zh = z(end); % z high

  warning('off','MATLAB:hg:UIControlSliderStepValueDifference');
  h_slider = uicontrol('Parent',h_f,'Style','slider','Position',[55,10,100,20],...
    'value',lsi, 'min',1, 'max',nl,'SliderStep',[1/max(nl-1,1) max(1/(nl-1),0.1)]);
  warning('on','MATLAB:hg:UIControlSliderStepValueDifference');
  h_centerbutton = uicontrol('Parent',h_f,'Style','pushbutton','String','Center','Position',[10,10,40,20]);
  string = [lambdatext ' = ' num2str(l(lsi),'%.3g')];
  h_slidertext = uicontrol('style','text','String',string,'BackgroundColor','w',...
    'horizontalAlignment','left','Position',[160,7,5.5*numel(string),20]);
  
  if numel(l) == 1 % Single-wavelength results
    set(h_slider,'Visible','off');
    set(h_centerbutton,'Visible','off');
    set(h_slidertext,'Visible','off');
    h_checkbox = uicontrol('Parent',h_f,'Style','checkbox','String','log plot','BackgroundColor','w','Position',[10,13,55,20]);
  else
    h_checkbox = uicontrol('Parent',h_f,'Style','checkbox','String','log plot','BackgroundColor','w','Position',[10,33,55,20]);
  end

  switch MCorFMC.boundaryType
    case 1
      NI_xneg_pad = MCorFMC.NI_xneg;
      NI_xneg_pad(ny+1,nz+1,1) = 0;
      NI_yneg_pad = MCorFMC.NI_yneg;
      NI_yneg_pad(nx+1,nz+1,1) = 0;
      NI_zneg_pad = MCorFMC.NI_zneg;
      NI_zneg_pad(nx+1,ny+1,1) = 0;
      NI_xpos_pad = MCorFMC.NI_xpos;
      NI_xpos_pad(ny+1,nz+1,1) = 0;
      NI_ypos_pad = MCorFMC.NI_ypos;
      NI_ypos_pad(nx+1,nz+1,1) = 0;
      NI_zpos_pad = MCorFMC.NI_zpos;
      NI_zpos_pad(nx+1,ny+1,1) = 0;
    
      h_f.UserData = {NI_xneg_pad, NI_yneg_pad, NI_zneg_pad, NI_xpos_pad, NI_ypos_pad, NI_zpos_pad};
    
      fprintf(['%.3g%% of ' fluorescenceOrIncident 'light hits the cuboid boundaries.\n'],100*(sum(MCorFMC.NI_xpos(:) + MCorFMC.NI_xneg(:))*dy*dz + sum(MCorFMC.NI_ypos(:) + MCorFMC.NI_yneg(:))*dx*dz + sum(MCorFMC.NI_zpos(:) + MCorFMC.NI_zneg(:))*dx*dy)/P_in);
      h_xneg = surface(repmat(xl,ny+1,nz+1),repmat(y',   1,nz+1),repmat(z ,ny+1,   1),NI_xneg_pad(:,:,lsi),'LineStyle','none');
      h_yneg = surface(repmat(x',   1,nz+1),repmat(yl,nx+1,nz+1),repmat(z ,nx+1,   1),NI_yneg_pad(:,:,lsi),'LineStyle','none');
      h_zneg = surface(repmat(x',   1,ny+1),repmat(y ,nx+1,   1),repmat(zl,nx+1,ny+1),NI_zneg_pad(:,:,lsi),'LineStyle','none');
      h_xpos = surface(repmat(xh,ny+1,nz+1),repmat(y',   1,nz+1),repmat(z ,ny+1,   1),NI_xpos_pad(:,:,lsi),'LineStyle','none');
      h_ypos = surface(repmat(x',   1,nz+1),repmat(yh,nx+1,nz+1),repmat(z ,nx+1,   1),NI_ypos_pad(:,:,lsi),'LineStyle','none');
      h_zpos = surface(repmat(x',   1,ny+1),repmat(y ,nx+1,   1),repmat(zh,nx+1,ny+1),NI_zpos_pad(:,:,lsi),'LineStyle','none');
      line(Lx/2*[1 1 1 1 -1 -1 -1 -1 1],Ly/2*[1 1 -1 -1 -1 -1 1 1 1],Lz*[1 0 0 1 1 0 0 1 1],'Color',[0.5 0.5 0.5]);
      line(Lx/2*[1 1 1 1 -1 -1 -1 -1 1],Ly/2*[1 1 -1 -1 -1 -1 1 1 1],Lz*[0 1 1 0 0 1 1 0 0],'Color',[0.5 0.5 0.5]);

      set(gca,'ZDir','reverse');
      xlabel('x [cm]');
      ylabel('y [cm]');
      zlabel('z [cm]');
      colormap(lincolormap);
      colorbar;
      axis tight
      axis equal
      set(gca,'fontsize',18)
      title(['Normalized ' fluorescenceOrNothing 'boundary irradiance [W/cm^2/W.incident]']);

      view(3)
      rotate3d on
      if ~verLessThan('matlab','9.0')
        setAxes3DPanAndZoomStyle(zoom(gca),gca,'camera');
      end
    case 2
      infwavecorrectionfactor = (size(MCorFMC.NI_zneg,1)/nx)^2;
      fprintf(['%.3g%% of ' fluorescenceOrIncident 'light hits the top cuboid boundary.\n'],100*sum(MCorFMC.NI_zneg(:))*dx*dy/P_in/infwavecorrectionfactor);

      h_f.UserData = {[],[],pagetranspose(MCorFMC.NI_zneg),[],[],[]};
      h_xneg = NaN;
      h_yneg = NaN;
      h_zneg = imagesc(size(MCorFMC.NI_zneg,1)/2*[-dx dx],size(MCorFMC.NI_zneg,2)/2*[-dy dy],MCorFMC.NI_zneg(:,:,lsi).');
      h_xpos = NaN;
      h_ypos = NaN;
      h_zpos = NaN;
      line(Lx/2*[-1 -1 1 1 -1],Ly/2*[-1 1 1 -1 -1],'Color',[1 1 1],'Linestyle','--');

      set(gca,'YDir','normal');
      xlabel('x [cm]');
      ylabel('y [cm]');
      colormap(lincolormap);
      colorbar;
      axis tight
      axis equal
      set(gca,'fontsize',18)
      title(['Normalized ' fluorescenceOrNothing 'boundary irradiance [W/cm^2/W.incident]']);
    case 3
      fprintf(['%.3g%% of ' fluorescenceOrIncident 'light hits the top or bottom cuboid boundaries.\n'],100*sum(MCorFMC.NI_zneg(:) + MCorFMC.NI_zpos(:))*dx*dy/P_in);

      h_f.UserData = {[],[],pagetranspose(MCorFMC.NI_zneg),[],[],pagetranspose(MCorFMC.NI_zpos)};
      h_xneg = NaN;
      h_yneg = NaN;
      subplot(2,1,1);
      h_zneg = imagesc(x,y,MCorFMC.NI_zneg(:,:,lsi).');
      h_xpos = NaN;
      h_ypos = NaN;
      set(gca,'YDir','normal');
      xlabel('x [cm]');
      ylabel('y [cm]');
      colormap(lincolormap);
      colorbar;
      axis tight
      axis equal
      set(gca,'fontsize',18)
      title(['Normalized ' fluorescenceOrNothing 'top boundary irradiance [W/cm^2/W.incident]']);

      subplot(2,1,2);
      h_zpos = imagesc(x,y,MCorFMC.NI_zpos(:,:,lsi).');

      set(gca,'YDir','normal');
      xlabel('x [cm]');
      ylabel('y [cm]');
      colormap(lincolormap);
      colorbar;
      axis tight
      axis equal
      set(gca,'fontsize',18)
      title(['Normalized ' fluorescenceOrNothing 'bottom boundary irradiance [W/cm^2/W.incident]']);
  end
  maxelement = -Inf;
  for iSurf = 1:numel(h_f.UserData)
    if ~isempty(h_f.UserData{iSurf})
      maxelement = max(maxelement,max(h_f.UserData{iSurf}(:)));
    end
  end
  if maxelement == 0
    maxelement = 1;
  end
  caxis([0 maxelement]);

  vars = struct('h_checkbox',h_checkbox,...
    'h_slider',h_slider,'h_centerbutton',h_centerbutton,...
    'h_slidertext',h_slidertext,'x',x,'y',y,'z',z,'l',l,...
    'h_xneg',h_xneg,'h_yneg',h_yneg,'h_zneg',h_zneg,...
    'h_xpos',h_xpos,'h_ypos',h_ypos,'h_zpos',h_zpos,...
    'linearColormap',lincolormap);

  set(h_checkbox,'Callback',{@redrawCuboidSurfaces,vars});
  set(h_slider,'Callback',{@redrawCuboidSurfaces,vars});
  set(h_centerbutton,'Callback',{@redrawCuboidSurfaces,vars});

  h_f.Name = ['Normalized ' fluorescenceOrNothing 'boundary irradiance'];
end