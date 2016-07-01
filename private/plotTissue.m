function plotTissue(T,tissueList,x,z,Nt,Nx,dx,Nz,dz)

imagesc(x,z,T,[1 Nt])
set(gca,'fontsize',18)
xlabel('x [cm]')
ylabel('z [cm]')
colorbar('YTickLabel',{'Escape','Air','Blood','Dermis', 'Epidermis',...
     'Skull','Grey matter','White matter','Hair', 'Gel'},'YTick',1:Nt);
cmap = makecmap(Nt);
colormap(cmap)
set(colorbar,'fontsize',1)
% label colorbar
xmin = min(x);
xmax = max(x);
zmin = min(z);
zmax = max(z);
zdiff = zmax-zmin;

for i=1:Nt
    yy = (Nt-i)/(Nt-1)*Nz*dz;
    text(Nx*dx*1.2,yy, tissueList(i).name,'fontsize',12)
end

text(xmax*0.9,zmin - zdiff*0.06, 'Tissue types','fontsize',18)
axis equal image
axis([xmin xmax zmin zmax])

return