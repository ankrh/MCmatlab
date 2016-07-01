function plotTemp(y,z,Temp_zy)

    image(y,z,Temp_zy)
    colorbar
    set(gca,'fontsize',18)
    xlabel('y [cm]')
    ylabel('z [cm]')
    colormap(makec2f)
    axis equal image

return