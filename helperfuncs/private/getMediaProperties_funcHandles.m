function MCinput = getMediaProperties_funcHandles(MCinput)
mediaProperties_raw = MCinput.G.mediaPropertiesFunc(MCinput.wavelength,MCinput.G.mediaPropParams); % The raw media properties, as defined by the user in the model file
mediaProperties_funcHandles = mediaProperties_raw; % The media properties, in which the char array expressions describing the dependences on intensity and temperature have been identified and converted to function handles

% Loop over all the media and throw error if the char arrays did not contain I or T
anyIntensityDependence = false;
anyTemperatureDependence = false;
for i=1:length(mediaProperties_funcHandles)
  if ischar(mediaProperties_funcHandles(i).mua)
    intensityDependence = contains(mediaProperties_funcHandles(i).mua,'I');
    anyIntensityDependence = anyIntensityDependence || intensityDependence;
    temperatureDependence = contains(mediaProperties_funcHandles(i).mua,'T');
    anyTemperatureDependence = anyTemperatureDependence || temperatureDependence;
    if ~intensityDependence && ~temperatureDependence
      error('Error: mua of %s depends on neither I (Intensity) nor T (Temperature)',mediaProperties_funcHandles(i).name);
    end
    mediaProperties_funcHandles(i).mua = str2func(['@(I,T)(' mediaProperties_funcHandles(i).mua ')']);
  end
  if ischar(mediaProperties_funcHandles(i).mus)
    intensityDependence = contains(mediaProperties_funcHandles(i).mus,'I');
    anyIntensityDependence = anyIntensityDependence || intensityDependence;
    temperatureDependence = contains(mediaProperties_funcHandles(i).mus,'T');
    anyTemperatureDependence = anyTemperatureDependence || temperatureDependence;
    if ~intensityDependence && ~temperatureDependence
      error('Error: mus of %s depends on neither I (Intensity) nor T (Temperature)',mediaProperties_funcHandles(i).name);
    end
    mediaProperties_funcHandles(i).mus = str2func(['@(I,T)(' mediaProperties_funcHandles(i).mus ')']);
  end
  if ischar(mediaProperties_funcHandles(i).g)
    intensityDependence = contains(mediaProperties_funcHandles(i).g,'I');
    anyIntensityDependence = anyIntensityDependence || intensityDependence;
    temperatureDependence = contains(mediaProperties_funcHandles(i).g,'T');
    anyTemperatureDependence = anyTemperatureDependence || temperatureDependence;
    if ~intensityDependence && ~temperatureDependence
      error('Error: g of %s depends on neither I (Intensity) nor T (Temperature)',mediaProperties_funcHandles(i).name);
    end
    mediaProperties_funcHandles(i).g = str2func(['@(I,T)(' mediaProperties_funcHandles(i).g ')']);
  end
end

MCinput.mediaProperties_funcHandles = mediaProperties_funcHandles;
MCinput.intensityDependent = anyIntensityDependence;
MCinput.temperatureDependent = anyTemperatureDependence;
end