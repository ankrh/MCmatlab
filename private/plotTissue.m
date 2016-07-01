function plotTissue(T,tissueList,x,z)

Nt = length(tissueList);
xmin = min(x);
xmax = max(x);
zmin = min(z);
zmax = max(z);
zdiff = zmax-zmin;

image(x,z,T)
xlabel('x [cm]')
ylabel('z [cm]')
axis equal image
cmap = makecmap(Nt);
colormap(cmap)
caxis([0 length(cmap)])
colorbar%('YTickLabel',{tissueList.name},'YTick',(1:Nt)-0.5);
set(gca,'fontsize',18)

text(xmax*0.9,zmin - zdiff*0.06, 'Tissue types','fontsize',18)

return