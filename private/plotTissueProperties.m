function plotTissueProperties(tissueList)

nT = length(tissueList);
cmap = colormap(hsv(nT));

subplot(231);
hold on;
for i=1:nT
    bar(i,tissueList(i).mua,'FaceColor',cmap(i,:));
end
set(gca,'XTick',1:nT,'XTickLabel',{tissueList.name},'XTickLabelRotation',45,'FontSize',12);
title('\mu_a [cm^{-1}]');
subplot(232);
hold on;
for i=1:nT
    bar(i,tissueList(i).mus,'FaceColor',cmap(i,:));
end
set(gca,'XTick',1:nT,'XTickLabel',{tissueList.name},'XTickLabelRotation',45,'FontSize',12);
title('\mu_s [cm^{-1}]');
subplot(233);
hold on;
for i=1:nT
    bar(i,tissueList(i).g,'FaceColor',cmap(i,:));
end
set(gca,'XTick',1:nT,'XTickLabel',{tissueList.name},'XTickLabelRotation',45,'FontSize',12);
title('g');
subplot(234);
hold on;
for i=1:nT
    bar(i,tissueList(i).VHC,'FaceColor',cmap(i,:));
end
set(gca,'XTick',1:nT,'XTickLabel',{tissueList.name},'XTickLabelRotation',45,'FontSize',12);
title('VHC [J/(cm^3*K)]');
subplot(235);
hold on;
for i=1:nT
    bar(i,tissueList(i).TC,'FaceColor',cmap(i,:));
end
set(gca,'XTick',1:nT,'XTickLabel',{tissueList.name},'XTickLabelRotation',45,'FontSize',12);
title('TC [W/(cm*K)]');

end