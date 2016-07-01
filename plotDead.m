function plotDead

directoryPath = './Data/';
load([directoryPath 'Input_spectrum'])

SC_thick = 0.002; %thickness of the stratum corneum layer [cm]
m = 0;
for n=F_min:dF:F_max
    m = m+1;
    myname = ['HeatSimOut_blood4_' num2str(m)];
    load([directoryPath myname '.mat']);

    p = 0;
    damage = 0;

    q = 0;
    bl_dam = 0;

    r = 0;
    wall_dam = 0;

    %% Plot showing the dead dermis and epidermis
    Tzy = squeeze(T(:,H_mci.Nx/2,:))'; % Make zy cut of the tissue matrix
    Temp_max_zy = squeeze(Temp_max(:,H_mci.Nx/2,:))';

    maxT = ones(length(Temp_max_zy(:,1,1)),length(Temp_max_zy(1,:,1)));
    Tdead = maxT;
    Temp_dead = 70;
    for i=1:length(Temp_max_zy(:,1,1))
        for ii = 1:length(Temp_max_zy(1,:,1))
            maxT(i,ii) = max(Temp_max_zy(i,ii,:));
            if (Tzy(i,ii) == 1 && maxT(i,ii)>Temp_dead) || (Tzy(i,ii) == 2 && maxT(i,ii)>Temp_dead) ...
                || (Tzy(i,ii) == 5 && maxT(i,ii)>Temp_dead)
                Tdead(i,ii) = 0;
            elseif Tzy(i,ii) == 3 && maxT(i,ii)>Temp_dead
                Tdead(i,ii) = 10;
            else
                Tdead(i,ii) = Tzy(i,ii);
            end

            if Tzy(i,ii) == 2 && maxT(i,ii) > Temp_dead
                p = p + 1;
                damage = damage + 1;
            elseif Tzy(i,ii) == 2 && maxT(i,ii) <= Temp_dead
                p = p + 1;
            end

            if Tzy(i,ii) == 3 && maxT(i,ii) > Temp_dead
                q = q + 1;
                bl_dam = bl_dam + 1;
            elseif Tzy(i,ii) == 3 && maxT(i,ii) <= Temp_dead
                q = q + 1;
            end

            if i>round((0.04-0.0020)/H_mci.dz) && i<=round((0.04)/H_mci.dz)
                if Tzy(i,ii) == 1 && maxT(i,ii) > Temp_dead
                    r = r + 1;
                    wall_dam = wall_dam + 1;
                elseif  Tzy(i,ii) == 1 && maxT(i,ii) <= Temp_dead
                    r = r + 1;
                end
            end
        end
    end

    dam_epi = damage/p*100; %ratio of epidermal damage [%]
    dam_blood = bl_dam/q*100; %ratio blood damage [%]
    dam_wall = wall_dam/r*100; %ratio of damage of the upper vessel wall [%]

    % At this point Tdead contains integers where 0=dead
    
    figure(1); clf
    plotTissue(Tzy,tissueList,x,z)

    fn = strrep(filename,'_',' ');
    title([num2str(fn) ' - ' num2str(n) ' J/cm2'])

    annotation('textbox', [0.01, 0.25, 0.2, 0.08], 'string',...
        ['Ratio of epidermal damage: ' num2str(dam_epi) ' %'],'FontSize',14)

    annotation('textbox', [0.01, 0.15, 0.2, 0.08], 'string',...
        ['Ratio of blood damage: ' num2str(dam_blood) ' %'],'FontSize',14)

    annotation('textbox', [0.01, 0.05, 0.2, 0.08], 'string',...
        ['Ratio of upper vessel wall damage: ' num2str(dam_wall) ' %'],'FontSize',14)

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
