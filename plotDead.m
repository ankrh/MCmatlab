function plotDead

directoryPath = './Data/';
load([directoryPath 'Input_spectrum'])

SC_thick = 0.002; %thickness of the stratum corneum layer [cm]
m = 0;
for n=F_min:dF:F_min
    m = m+1;
    num = num2str(m);
    numb = num2str(n);
    number = strrep(numb,'.',',');
    load([directoryPath 'HeatSimOut_blood4_' num '.mat']);
%     load(['VL2_bund_' number '.mat']);

Tzy = reshape(T(:,Nx/2,:),Ny,Nz)'; % Make zy cut of the tissue matrix

%% Tissue model plot
% At this point Tzy contains integers where 4=dermis,
% 5=epidermis, 9=blood and 10=gel.
% This is changed to, 1=dermis, 2=epidermis, 3=blood, 4=gel, 5=stratum corneum.
Tzy(Tzy==4)=1; Tzy(Tzy==5)=2; Tzy(Tzy==9)=3; Tzy(Tzy==10)=4;

for iz = 1:size(Tzy,1)
    if iz>round((0.02)/dz) && iz<=round((0.02+SC_thick)/dz)
        Tzy(iz,:) = 5;
    end
end

% Make the color map
Nt = 5; % number of tissues (including dead dermis)
cmap = zeros(64,3);
dj = 0.05;
for i=1:64
    j = round((i-dj)/64*(Nt-1));
    if  j<=1-dj, cmap(i,:) = [1 0.8 0.8]; % dermis
    elseif  j<=2-dj, cmap(i,:) = [0.5 0.2 0.2]; % epidermis
    elseif  j<=3-dj, cmap(i,:) = [1 0 0]; % blood
    elseif  j<=4-dj, cmap(i,:) = [0 0 1]; % gel
    elseif  j<=5-dj, cmap(i,:) = [0.8 0.7 0.8]; %SC
    end
end

% figure('Position',[10 50 800 600]);clf
% imagesc(y,z,Tzy,[1 Nt])
% hold on
% colormap(cmap)
% colorbar('YTickLabel',{'Dermis', 'Epidermis',...
%     'Blood', 'Gel','Stratum Corneum'},'YTick',1:Nt);
% set(gca,'fontsize',18)
% xlabel('y [cm]')
% ylabel('z [cm]')
% title('Tissue Model')
% axis equal image

   
p = 0;
damage = 0;

q = 0;
bl_dam = 0;

r = 0;
wall_dam = 0;

%% Plot showing the dead dermis and epidermis
maxT = ones(length(Tempzy(:,1,1)),length(Tempzy(1,:,1)));
Tdead = maxT;
Temp_dead = 70;
for i=1:length(Tempzy(:,1,1))
    for ii = 1:length(Tempzy(1,:,1))
        maxT(i,ii) = max(Tempzy(i,ii,:));
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
                
        if i>round((0.04-0.0020)/dz) && i<=round((0.04)/dz)
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

% At this point Tdead contains integers where 0=dead, 4=dermis,
% 5=epidermis, 9=hair and 10=gel.
% This is changed to 0=dead, 1=dermis, 2=epidermis, 3=hair and 4=gel.
% Tdead(Tdead==4)=1; Tdead(Tdead==5)=2; Tdead(Tdead==9)=3; Tdead(Tdead==10)=4;

% Make the color map
Nt = 7; % number of tissues (including dead dermis)
cmap = zeros(64,3);
dj = 0.05;
for i=1:64
    j = round((i-dj)/64*(Nt-1));
    if      j<=1-dj, cmap(i,:) = [1 1 0]; % dead
    elseif  j<=2-dj, cmap(i,:) = [1 0.8 0.8]; % dermis
    elseif  j<=3-dj, cmap(i,:) = [0.5 0.2 0.2]; % epidermis
    elseif  j<=4-dj, cmap(i,:) = [1 0 0]; % Blood
    elseif  j<=5-dj, cmap(i,:) = [0 0 1]; % gel
    elseif  j<=6-dj, cmap(i,:) = [0.8 0.7 0.8]; %SC
    elseif  j<=10-dj, cmap(i,:) = [0,0,0]; %dead blood
    end
end

fn = strrep(filename,'_',' ');

figure('units','normalized','outerposition',[0 0 0.95 1]); clf
imagesc(y,z,Tdead,[0 Nt-1])
hold on
colormap(cmap)
colorbar('YTickLabel',{'Dead', 'Dermis', 'Epidermis',...
    'Blood', 'Gel','Stratum Corneum','DeadBlood'},'YTick',0:Nt-1);
set(gca,'fontsize',18)
xlabel('y [cm]')
ylabel('z [cm]')
title([num2str(fn) ' - ' num2str(n) ' J/cm2'])
axis equal image

annotation('textbox', [0.01, 0.25, 0.2, 0.08], 'string',...
    ['Ratio of epidermal damage: ' num2str(dam_epi) ' %'],'FontSize',14)

annotation('textbox', [0.01, 0.15, 0.2, 0.08], 'string',...
    ['Ratio of blood damage: ' num2str(dam_blood) ' %'],'FontSize',14)

annotation('textbox', [0.01, 0.05, 0.2, 0.08], 'string',...
    ['Ratio of upper vessel wall damage: ' num2str(dam_wall) ' %'],'FontSize',14)


% yint = 45:55; zint = 50:70;
% figure('Position',[10 50 600 600]);clf
% imagesc(y(yint),z(zint),Tdead(zint,yint),[0 Nt-1])
% hold on
% colormap(cmap)
% colorbar('YTickLabel',{'Dead', 'Dermis', 'Epidermis',...
%     'Blood', 'Gel','DeadBlood'},'YTick',0:Nt-1);
% set(gca,'fontsize',18)
% xlabel('y [cm]')
% ylabel('z [cm]')
% title('Dead dermis/epidermis')
% axis equal image

end
return
%% Animate Temperature
clear TempMovie
t_int = length(Tempzy(1,1,:));
time = 1:t_int; % in milliseconds

% vid = VideoWriter('Temperature_Movie_630nm.avi');
% vid.Quality = 100;
% vid.FrameRate = 2;
% open(vid);
h = figure;
set(h, 'Renderer', 'painters')

for i=1:45%round(t_int/2)
    imagesc(y,z,Tempzy(:,:,i))
    hold on
    text(max(x)*0.9,min(z)-0.04*max(z),'T [^{\circ}C]','fontsize',18)
    colorbar
    set(gca,'fontsize',18)
    xlabel('y [cm]')
    ylabel('z [cm]')
    title(['Time = ' num2str(time(i)) ' ms'])
    axis equal image
    drawnow
    %writeVideo(vid, getframe(h));
    TempMovie(i) = getframe(h);
end
%close(vid);
%winopen('Temperature_Movie_630nm.avi')
    movie2avi(TempMovie, 'Temperature_Movie_630nm','compression','None','fps',2)
return
%% Old code


figure;clf
imagesc(y,z,maxT)
hold on
text(max(x)*0.9,min(z)-0.04*max(z),'T [^{\circ}C]','fontsize',18)
colorbar
set(gca,'fontsize',18)
xlabel('y [cm]')
ylabel('z [cm]')
title('Temperature  [^{\circ}C] ')
%colormap(makec2f)
axis equal image

maxTemp_post_light = max(max(max(Temp_post_light)))
[maxTemp_papilla,I] = max(Tempzy(63,50,:)) %maximum temperature in the center of the papilla

Temp_post_lightzy = reshape(Temp_post_light(:,Nx/2,:),Ny,Nz)';


figure(1);clf
imagesc(y,z,Tzy,[1 Nt])
hold on
cmap = makecmap(Nt);
colormap(cmap)
colorbar
set(gca,'fontsize',18)
set(colorbar,'fontsize',1)
xlabel('x [cm]')
ylabel('z [cm]')
title('Tissue types')
for i=1:Nt
    yy = min(z) + (Nt-i)/(Nt-1)*max(z);
    text(max(x)+1/6*(max(x)-min(x)),yy, sprintf('%d %s',i,tissue(i).s),'fontsize',12)
end
axis equal image

figure(5);clf
imagesc(y,z,Temp_post_lightzy)
hold on
text(max(x)*0.9,min(z)-0.04*max(z),'T [^{\circ}C]','fontsize',18)
colorbar
set(gca,'fontsize',18)
xlabel('y [cm]')
ylabel('z [cm]')
title('Temperature  [^{\circ}C] ')
%colormap(makec2f)
axis equal image

figure(6);clf
imagesc(y,z,Tempzy(:,:,30))
hold on
text(max(x)*0.9,min(z)-0.04*max(z),'T [^{\circ}C]','fontsize',18)
colorbar
set(gca,'fontsize',18)
xlabel('y [cm]')
ylabel('z [cm]')
title('Temperature  [^{\circ}C] ')
%colormap(makec2f)
axis equal image

figure(7);clf
imagesc(y,z,Temp_post_lightzy(:,:)-Tempzy(:,:,1))
hold on
text(max(x)*0.9,min(z)-0.04*max(z),'T [^{\circ}C]','fontsize',18)
colorbar
set(gca,'fontsize',18)
xlabel('y [cm]')
ylabel('z [cm]')
title('Temperature  [^{\circ}C] ')
%colormap(makec2f)
axis equal image
%print -djpeg -r300 'Fig_Azy.jpg'