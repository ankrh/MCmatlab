function plotTissue(T,tissueList,x,z)

Nt = length(tissueList);
Nx = length(x);
xmin = min(x);
xmax = max(x);
Nz = length(z);
zmin = min(z);
zmax = max(z);
zdiff = zmax-zmin;

imagesc(x,z,T,[1 Nt])
xlabel('x [cm]')
ylabel('z [cm]')
axis([xmin xmax zmin zmax])
axis equal image
cmap = makecmap(Nt);
colormap(cmap)
colorbar('YTickLabel',{tissueList.name},'YTick',1:Nt);
set(gca,'fontsize',18)

text(xmax*0.9,zmin - zdiff*0.06, 'Tissue types','fontsize',18)

return