function plotMediaProperties(mediaProperties)

clf;

nM = length(mediaProperties);
cmap = colormap(lines(nM));

plotThermalProperties = (isfield(mediaProperties,'VHC') && isfield(mediaProperties,'TC') && length([mediaProperties.VHC]) == nM && length([mediaProperties.TC]) == nM);
plotFluorescenceProperties = any([mediaProperties.Y] > 0);

if plotThermalProperties && plotFluorescenceProperties
    rows = 3;
elseif plotThermalProperties || plotFluorescenceProperties
    rows = 2;
else
    rows = 1;
end

subplot(rows,3,1);
hold on;
for i=1:nM
    bar(i,mediaProperties(i).mua,'FaceColor',cmap(i,:));
end
set(gca,'XTick',1:nM,'XTickLabel',{mediaProperties.name},'XTickLabelRotation',45,'FontSize',12,'Box','on','YGrid','on','YMinorGrid','on');
title('\mu_a [cm^{-1}]');

subplot(rows,3,2);
hold on;
for i=1:nM
    bar(i,mediaProperties(i).mus,'FaceColor',cmap(i,:));
end
set(gca,'XTick',1:nM,'XTickLabel',{mediaProperties.name},'XTickLabelRotation',45,'FontSize',12,'Box','on','YGrid','on','YMinorGrid','on');
title('\mu_s [cm^{-1}]');

subplot(rows,3,3);
hold on;
for i=1:nM
    bar(i,mediaProperties(i).g,'FaceColor',cmap(i,:));
end
set(gca,'XTick',1:nM,'XTickLabel',{mediaProperties.name},'XTickLabelRotation',45,'FontSize',12,'Box','on','YGrid','on','YMinorGrid','on');
title('g');

if plotThermalProperties
    subplot(rows,3,4);
    hold on;
    for i=1:nM
        bar(i,mediaProperties(i).VHC,'FaceColor',cmap(i,:));
    end
    set(gca,'XTick',1:nM,'XTickLabel',{mediaProperties.name},'XTickLabelRotation',45,'FontSize',12,'Box','on','YGrid','on','YMinorGrid','on');
    title('VHC [J/(cm^3*K)]');

    subplot(rows,3,5);
    hold on;
    for i=1:nM
        bar(i,mediaProperties(i).TC,'FaceColor',cmap(i,:));
    end
    set(gca,'XTick',1:nM,'XTickLabel',{mediaProperties.name},'XTickLabelRotation',45,'FontSize',12,'Box','on','YGrid','on','YMinorGrid','on');
    title('TC [W/(cm*K)]');

    if  plotFluorescenceProperties
        subplot(rows,3,7);
        hold on;
        for i=1:nM
            bar(i,mediaProperties(i).Y,'FaceColor',cmap(i,:));
        end
        set(gca,'XTick',1:nM,'XTickLabel',{mediaProperties.name},'XTickLabelRotation',45,'FontSize',12,'Box','on','YGrid','on','YMinorGrid','on');
        title('\epsilon');

        subplot(rows,3,8);
        hold on;
        for i=1:nM
            bar(i,mediaProperties(i).sat,'FaceColor',cmap(i,:));
        end
        set(gca,'XTick',1:nM,'XTickLabel',{mediaProperties.name},'XTickLabelRotation',45,'FontSize',12,'Box','on','YGrid','on','YMinorGrid','on');
        title('Fluor. sat. [W/cm^3]');
    end
elseif  plotFluorescenceProperties
    subplot(rows,3,4);
    hold on;
    for i=1:nM
        bar(i,mediaProperties(i).Y,'FaceColor',cmap(i,:));
    end
    set(gca,'XTick',1:nM,'XTickLabel',{mediaProperties.name},'XTickLabelRotation',45,'FontSize',12,'Box','on','YGrid','on','YMinorGrid','on');
    title('Y');

    subplot(rows,3,5);
    hold on;
    for i=1:nM
        bar(i,mediaProperties(i).sat,'FaceColor',cmap(i,:));
    end
    set(gca,'XTick',1:nM,'XTickLabel',{mediaProperties.name},'XTickLabelRotation',45,'FontSize',12,'Box','on','YGrid','on','YMinorGrid','on');
    title('Fluor. sat. [W/cm^2]');
end
end