function h_f = plotMediaProperties(nFig,model,simType)

switch simType
  case 1
    mP = model.MC.mediaProperties;
    mP_fH = model.MC.mediaProperties_funcHandles;
  case 2
    mP = model.FMC.mediaProperties;
    mP_fH = model.FMC.mediaProperties_funcHandles;
  case 3
    mP = model.HS.mediaProperties;
    mP_fH = model.HS.mediaProperties_funcHandles;
end

if(~ishandle(nFig))
  h_f = figure(nFig);
  h_f.Position = [40 160 1100 650];
else
  h_f = figure(nFig);
end
h_f.Color = 'w';

clf;

nM = length(mP_fH);
cmap = colormap(lines);

if simType <= 2
  subplot(2,2,1);
  hold on;
  barIdx = 1;
  splitMediumIdx = 1;
  nameCellArray = {};
  for mediumIdx=1:nM
    if isa(mP_fH(mediumIdx).mua,'function_handle')
      bar(barIdx,min([mP(splitMediumIdx:(splitMediumIdx + mP_fH(mediumIdx).nBins - 1)).mua]),'FaceColor',cmap(mediumIdx,:));
      nameCellArray(barIdx) = {[mP_fH(mediumIdx).name ' (min)']};
      barIdx = barIdx+1;
      bar(barIdx,max([mP(splitMediumIdx:(splitMediumIdx + mP_fH(mediumIdx).nBins - 1)).mua]),'FaceColor',cmap(mediumIdx,:));
      nameCellArray(barIdx) = {[mP_fH(mediumIdx).name ' (max)']};
      barIdx = barIdx+1;
    else
      bar(barIdx,mP(splitMediumIdx).mua,'FaceColor',cmap(barIdx,:));
      nameCellArray(barIdx) = {mP_fH(mediumIdx).name};
      barIdx = barIdx+1;
    end
    if ~isa(mP_fH(mediumIdx).mua,'function_handle')
      splitMediumIdx = splitMediumIdx + 1;
    else
      splitMediumIdx = splitMediumIdx + mP_fH(mediumIdx).nBins;
    end
  end
  set(gca,'XTick',1:(barIdx-1),'XTickLabel',nameCellArray,'XTickLabelRotation',45,'FontSize',12,'Box','on','YGrid','on','YMinorGrid','on');
  title('\mu_a [cm^{-1}]');

  subplot(2,2,2);
  hold on;
  barIdx = 1;
  splitMediumIdx = 1;
  nameCellArray = {};
  for mediumIdx=1:nM
    if isa(mP_fH(mediumIdx).mus,'function_handle')
      bar(barIdx,min([mP(splitMediumIdx:(splitMediumIdx + mP_fH(mediumIdx).nBins - 1)).mus]),'FaceColor',cmap(mediumIdx,:));
      nameCellArray(barIdx) = {[mP_fH(mediumIdx).name ' (min)']};
      barIdx = barIdx+1;
      bar(barIdx,max([mP(splitMediumIdx:(splitMediumIdx + mP_fH(mediumIdx).nBins - 1)).mus]),'FaceColor',cmap(mediumIdx,:));
      nameCellArray(barIdx) = {[mP_fH(mediumIdx).name ' (max)']};
      barIdx = barIdx+1;
    else
      bar(barIdx,mP(splitMediumIdx).mus,'FaceColor',cmap(barIdx,:));
      nameCellArray(barIdx) = {mP_fH(mediumIdx).name};
      barIdx = barIdx+1;
    end
    if ~isa(mP_fH(mediumIdx).mus,'function_handle')
      splitMediumIdx = splitMediumIdx + 1;
    else
      splitMediumIdx = splitMediumIdx + mP_fH(mediumIdx).nBins;
    end
  end
  set(gca,'XTick',1:(barIdx-1),'XTickLabel',nameCellArray,'XTickLabelRotation',45,'FontSize',12,'Box','on','YGrid','on','YMinorGrid','on');
  title('\mu_s [cm^{-1}]');

  subplot(2,2,3);
  hold on;
  barIdx = 1;
  splitMediumIdx = 1;
  nameCellArray = {};
  for mediumIdx=1:nM
    if isa(mP_fH(mediumIdx).g,'function_handle')
      bar(barIdx,min([mP(splitMediumIdx:(splitMediumIdx + mP_fH(mediumIdx).nBins - 1)).g]),'FaceColor',cmap(mediumIdx,:));
      nameCellArray(barIdx) = {[mP_fH(mediumIdx).name ' (min)']};
      barIdx = barIdx+1;
      bar(barIdx,max([mP(splitMediumIdx:(splitMediumIdx + mP_fH(mediumIdx).nBins - 1)).g]),'FaceColor',cmap(mediumIdx,:));
      nameCellArray(barIdx) = {[mP_fH(mediumIdx).name ' (max)']};
      barIdx = barIdx+1;
    else
      bar(barIdx,mP(splitMediumIdx).g,'FaceColor',cmap(barIdx,:));
      nameCellArray(barIdx) = {mP_fH(mediumIdx).name};
      barIdx = barIdx+1;
    end
    if ~isa(mP_fH(mediumIdx).g,'function_handle')
      splitMediumIdx = splitMediumIdx + 1;
    else
      splitMediumIdx = splitMediumIdx + mP_fH(mediumIdx).nBins;
    end
  end
  set(gca,'XTick',1:(barIdx-1),'XTickLabel',nameCellArray,'XTickLabelRotation',45,'FontSize',12,'Box','on','YGrid','on','YMinorGrid','on');
  title('g');

  subplot(2,2,4);
  hold on;
  if (simType == 1 && ~model.MC.matchedInterfaces) || ...
     (simType == 2 && ~model.FMC.matchedInterfaces)
    for i=1:nM
      bar(i,mP_fH(i).n,'FaceColor',cmap(i,:));
    end
    set(gca,'XTick',1:nM,'XTickLabel',{mP_fH.name},'XTickLabelRotation',45,'FontSize',12,'Box','on','YGrid','on','YMinorGrid','on');
  else
    ax = gca;
    ax.XAxis.Visible = 'off';
    ax.YAxis.Visible = 'off';
    str = {'Assuming index','matched interfaces:','n = 1'};
    text(0.5,0.5,str,'FontSize',12,'HorizontalAlignment','center');
  end
  title('n','FontSize',12);
else
  subplot(1,2,1);
  hold on;
  barIdx = 1;
  splitMediumIdx = 1;
  nameCellArray = {};
  for mediumIdx=1:nM
    if isa(mP_fH(mediumIdx).VHC,'function_handle')
      bar(barIdx,min([mP(splitMediumIdx:(splitMediumIdx + mP_fH(mediumIdx).nBins - 1)).VHC]),'FaceColor',cmap(mediumIdx,:));
      nameCellArray(barIdx) = {[mP_fH(mediumIdx).name ' (min)']};
      barIdx = barIdx+1;
      bar(barIdx,max([mP(splitMediumIdx:(splitMediumIdx + mP_fH(mediumIdx).nBins - 1)).VHC]),'FaceColor',cmap(mediumIdx,:));
      nameCellArray(barIdx) = {[mP_fH(mediumIdx).name ' (max)']};
      barIdx = barIdx+1;
    else
      bar(barIdx,mP(splitMediumIdx).VHC,'FaceColor',cmap(barIdx,:));
      nameCellArray(barIdx) = {mP_fH(mediumIdx).name};
      barIdx = barIdx+1;
    end
    if ~isa(mP_fH(mediumIdx).VHC,'function_handle')
      splitMediumIdx = splitMediumIdx + 1;
    else
      splitMediumIdx = splitMediumIdx + mP_fH(mediumIdx).nBins;
    end
  end
  set(gca,'XTick',1:(barIdx-1),'XTickLabel',nameCellArray,'XTickLabelRotation',45,'FontSize',12,'Box','on','YGrid','on','YMinorGrid','on');
  title('Volumetric Heat Capacity [J/(cm^3*K)]');

  subplot(1,2,2);
  hold on;
  barIdx = 1;
  splitMediumIdx = 1;
  nameCellArray = {};
  for mediumIdx=1:nM
    if isa(mP_fH(mediumIdx).TC,'function_handle')
      bar(barIdx,min([mP(splitMediumIdx:(splitMediumIdx + mP_fH(mediumIdx).nBins - 1)).TC]),'FaceColor',cmap(mediumIdx,:));
      nameCellArray(barIdx) = {[mP_fH(mediumIdx).name ' (min)']};
      barIdx = barIdx+1;
      bar(barIdx,max([mP(splitMediumIdx:(splitMediumIdx + mP_fH(mediumIdx).nBins - 1)).TC]),'FaceColor',cmap(mediumIdx,:));
      nameCellArray(barIdx) = {[mP_fH(mediumIdx).name ' (max)']};
      barIdx = barIdx+1;
    else
      bar(barIdx,mP(splitMediumIdx).TC,'FaceColor',cmap(barIdx,:));
      nameCellArray(barIdx) = {mP_fH(mediumIdx).name};
      barIdx = barIdx+1;
    end
    if ~isa(mP_fH(mediumIdx).TC,'function_handle')
      splitMediumIdx = splitMediumIdx + 1;
    else
      splitMediumIdx = splitMediumIdx + mP_fH(mediumIdx).nBins;
    end
  end
  set(gca,'XTick',1:(barIdx-1),'XTickLabel',nameCellArray,'XTickLabelRotation',45,'FontSize',12,'Box','on','YGrid','on','YMinorGrid','on');
  title('Thermal Conductivity [W/(cm*K)]');
end
end