function plotTissue(T,tissueList,x,z)

Nt = length(tissueList);
Nx = length(x);
xmin = min(x);
xmax = max(x);
dx = (xmax-xmin)/Nx;
Nz = length(z);
zmin = min(z);
zmax = max(z);
zdiff = zmax-zmin;
dz = zdiff/Nz;

imagesc(x,z,T,[1 Nt])
set(gca,'fontsize',18)
xlabel('x [cm]')
ylabel('z [cm]')
colorbar('YTickLabel',{tissueList.name},'YTick',1:Nt);
cmap = makecmap(Nt);
colormap(cmap)
set(colorbar,'fontsize',1)

for i=1:Nt
    yy = (Nt-i)/(Nt-1)*Nz*dz;
    text(Nx*dx*1.2,yy, tissueList(i).name,'fontsize',12)
end

text(xmax*0.9,zmin - zdiff*0.06, 'Tissue types','fontsize',18)
axis equal image
axis([xmin xmax zmin zmax])

return