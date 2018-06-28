function plotTissueProperties(tissueList)

clf;

nT = length(tissueList);
cmap = colormap(lines(nT));

plotThermalProperties = (isfield(tissueList,'VHC') && isfield(tissueList,'TC') && length([tissueList.VHC]) == nT && length([tissueList.TC]) == nT);
plotFluorescenceProperties = any([tissueList.Y] > 0);

if plotThermalProperties && plotFluorescenceProperties
    rows = 3;
elseif plotThermalProperties || plotFluorescenceProperties
    rows = 2;
else
    rows = 1;
end

subplot(rows,3,1);
hold on;
for i=1:nT
    bar(i,tissueList(i).mua,'FaceColor',cmap(i,:));
end
set(gca,'XTick',1:nT,'XTickLabel',{tissueList.name},'XTickLabelRotation',45,'FontSize',12);
title('\mu_a [cm^{-1}]');

subplot(rows,3,2);
hold on;
for i=1:nT
    bar(i,tissueList(i).mus,'FaceColor',cmap(i,:));
end
set(gca,'XTick',1:nT,'XTickLabel',{tissueList.name},'XTickLabelRotation',45,'FontSize',12);
title('\mu_s [cm^{-1}]');

subplot(rows,3,3);
hold on;
for i=1:nT
    bar(i,tissueList(i).g,'FaceColor',cmap(i,:));
end
set(gca,'XTick',1:nT,'XTickLabel',{tissueList.name},'XTickLabelRotation',45,'FontSize',12);
title('g');

if plotThermalProperties
    subplot(rows,3,4);
    hold on;
    for i=1:nT
        bar(i,tissueList(i).VHC,'FaceColor',cmap(i,:));
    end
    set(gca,'XTick',1:nT,'XTickLabel',{tissueList.name},'XTickLabelRotation',45,'FontSize',12);
    title('VHC [J/(cm^3*K)]');
    
    subplot(rows,3,5);
    hold on;
    for i=1:nT
        bar(i,tissueList(i).TC,'FaceColor',cmap(i,:));
    end
    set(gca,'XTick',1:nT,'XTickLabel',{tissueList.name},'XTickLabelRotation',45,'FontSize',12);
    title('TC [W/(cm*K)]');
    
    if  plotFluorescenceProperties
        subplot(rows,3,7);
        hold on;
        for i=1:nT
            bar(i,tissueList(i).Y,'FaceColor',cmap(i,:));
        end
        set(gca,'XTick',1:nT,'XTickLabel',{tissueList.name},'XTickLabelRotation',45,'FontSize',12);
        title('\epsilon');
        
        subplot(rows,3,8);
        hold on;
        for i=1:nT
            bar(i,tissueList(i).sat,'FaceColor',cmap(i,:));
        end
        set(gca,'XTick',1:nT,'XTickLabel',{tissueList.name},'XTickLabelRotation',45,'FontSize',12);
        title('Fluor. sat. [W/cm^3]');
    end
elseif  plotFluorescenceProperties
    subplot(rows,3,4);
    hold on;
    for i=1:nT
        bar(i,tissueList(i).Y,'FaceColor',cmap(i,:));
    end
    set(gca,'XTick',1:nT,'XTickLabel',{tissueList.name},'XTickLabelRotation',45,'FontSize',12);
    title('Y');

    subplot(rows,3,5);
    hold on;
    for i=1:nT
        bar(i,tissueList(i).sat,'FaceColor',cmap(i,:));
    end
    set(gca,'XTick',1:nT,'XTickLabel',{tissueList.name},'XTickLabelRotation',45,'FontSize',12);
    title('Fluor. sat. [W/cm^2]');
end
end