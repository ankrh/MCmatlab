function plotDead

directoryPath = './Data/';
load([directoryPath 'Input_spectrum'])

m = 0;
for n=F_min:dF:F_max
    m = m+1;
    myname = ['HeatSimOut_blood4_' num2str(m)];
    load([directoryPath myname '.mat']);
    
    Temp_dead = 70;

    Tzy = squeeze(T(:,H_mci.Nx/2,:))'; % Make zy cut of the tissue matrix
    Temp_max_zy = squeeze(Temp_max(:,H_mci.Nx/2,:))';

    tissueIndex_epidermis = find(strcmp({tissueList.name},'epidermis'),1);
    total_epidermis = sum(Tzy==tissueIndex_epidermis);
    damage_epidermis = sum(Tzy==tissueIndex_epidermis & Temp_max_zy >= Temp_dead);

    tissueIndex_blood = find(strcmp({tissueList.name},'blood'),1);
    total_blood = sum(Tzy==tissueIndex_blood);
    damage_blood = sum(Tzy==tissueIndex_blood & Temp_max_zy >= Temp_dead);

    dam_epi = damage_epidermis/total_epidermis*100; %ratio of epidermal damage [%]
    dam_blood = damage_blood/total_blood*100; %ratio blood damage [%]

    %% Plot showing the dead dermis and epidermis

    Tzy(Temp_max_zy >= Temp_dead) = 0;
    % At this point Tdead contains integers where 0=dead
    
    figure(1); clf
    plotTissue(Tzy,tissueList,x,z)

    fn = strrep(filename,'_',' ');
    title([num2str(fn) ' - ' num2str(n) ' J/cm2'])

    annotation('textbox', [0.01, 0.25, 0.2, 0.08], 'string',...
        ['Ratio of epidermal damage: ' num2str(dam_epi) ' %'],'FontSize',14)

    annotation('textbox', [0.01, 0.15, 0.2, 0.08], 'string',...
        ['Ratio of blood damage: ' num2str(dam_blood) ' %'],'FontSize',14)

    Temp_post_light_zy = squeeze(Temp_post_light(:,H_mci.Nx/2,:))';
    figure(5);clf
    image(y,z,Temp_post_light_zy)
    hold on
    text(max(x)*0.9,min(z)-0.04*max(z),'T [^{\circ}C]','fontsize',18)
    colorbar
    set(gca,'fontsize',18)
    xlabel('y [cm]')
    ylabel('z [cm]')
    title('Temperature after Illumination [^{\circ}C] ')
    colormap(makec2f)
    axis equal image
    name = sprintf('%s%s_T_post_light_zy.jpg',directoryPath,myname);
    print('-djpeg','-r300',name)

    Temp_post_diffuse_zy = squeeze(Temp_post_diffuse(:,H_mci.Nx/2,:))';
    figure(6);clf
    image(y,z,Temp_post_diffuse_zy)
    hold on
    text(max(x)*0.9,min(z)-0.04*max(z),'T [^{\circ}C]','fontsize',18)
    colorbar
    set(gca,'fontsize',18)
    xlabel('y [cm]')
    ylabel('z [cm]')
    title('Temperature after Diffusion [^{\circ}C]')
    colormap(makec2f)
    axis equal image
    name = sprintf('%s%s_T_post_diffuse_zy.jpg',directoryPath,myname);
    print('-djpeg','-r300',name)

    figure(7);clf
    image(y,z,Temp_max_zy)
    hold on
    text(max(x)*0.9,min(z)-0.04*max(z),'T [^{\circ}C]','fontsize',18)
    colorbar
    set(gca,'fontsize',18)
    xlabel('y [cm]')
    ylabel('z [cm]')
    title('maximum Temperature reached [^{\circ}C]')
    colormap(makec2f)
    axis equal image
    name = sprintf('%s%s_T_max_zy.jpg',directoryPath,myname);
    print('-djpeg','-r300',name)

end

return
