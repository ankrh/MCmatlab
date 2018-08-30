function plotMCmatlabHeat(HSinput,HSoutput)
%   Created 2018 by Dominik Marti and Anders K. Hansen, DTU Fotonik
%
%   Requires
%       plotVolumetric.m
%

G = HSinput.G;
nM = length(G.mediaProperties); % Number of different media in simulation
numTemperatureSensors = size(HSinput.tempSensorPositions,1);

%% Plot the geometry showing the temperature sensor locations and the sensor data
if(numTemperatureSensors)
    geometryFigure = plotVolumetric(22,G.x,G.y,G.z,G.M,'MCmatlab_GeometryIllustration',G.mediaProperties,'slicePositions',HSinput.slicePositions);
    geometryFigure.Name = 'Temperature sensor illustration';
    title('Temperature sensor illustration');

    for i=numTemperatureSensors:-1:1
        sensorLabels{i,1} = num2str(i);
    end
    text(HSinput.tempSensorPositions(:,1),HSinput.tempSensorPositions(:,2),HSinput.tempSensorPositions(:,3),sensorLabels,'HorizontalAlignment','center','FontSize',18);

    if(~ishandle(23))
        temperatureSensorFigure = figure(23);
        temperatureSensorFigure.Position = [40 80 1100 650];
    else
        temperatureSensorFigure = figure(23);
    end
    clf;
    temperatureSensorFigure.Name = 'Temperature sensors';
    plot(HSoutput.sensorsTimeVector,HSoutput.sensorTemps,'LineWidth',2);
    set(gca,'FontSize',16);
    xlabel('Time [sec]');
    ylabel('Temperature [deg C]');
    title('Temperature sensors');
    xlim(HSoutput.sensorsTimeVector([1 end]));
    legend(sensorLabels,'Location','best');
    grid on;grid minor;
end

%% Plot thermal damage
if ~isnan(HSoutput.Omega(1))
    M_damage = G.M;
    M_damage(HSoutput.Omega > 1) = nM + 1;
    G.mediaProperties(nM + 1).name = 'damage';
    damageFigure = plotVolumetric(25,G.x,G.y,G.z,M_damage,'MCmatlab_GeometryIllustration',G.mediaProperties,'slicePositions',HSinput.slicePositions);
    damageFigure.Name = 'Thermal damage illustration';
    title('Thermal damage illustration');
    fprintf('%.2e cm^3 was thermally damaged.\n',G.dx*G.dy*G.dz*sum(sum(sum(HSoutput.Omega > 1))));
end
drawnow;
end

