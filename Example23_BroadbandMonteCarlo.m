%% Decription
% This example showcases the ability to launch broadband simulations, in
% which multiple wavelengths of excitation and fluorescence light are
% simulated. The geometry is loosely based on that of example 4 (blood
% vessel), but with a sphere of some fluorescing stained cancerous tissue
% inserted into the dermis and a droplet of water with Lucifer Yellow CH on
% top of the epidermis. The use of the wavelength-dependent calc_mua() and
% calc_mus() functions in the media properties allows us to easily take
% into account the tissues' absorption and scattering variations with
% wavelength.
% 
% Note the lambda slider that is now present in the output volumetric plots.
% 
% Each excitation wavelength has a power determined by the user-defined
% function handle stored in the model.MC.spectrumFunc property. This
% function can be defined along with the geometry function and the media
% properties function in the model file itself as shown below and must take
% one input argument (the wavelength) and return one output argument (the
% spectral power). If no function handle is assigned to spectrumFunc, it
% will be assumed that all wavelengths share the power equally.
% 
% The spectrumFunc values do not need to be normalized, as MCmatlab
% normalizes the powers at all the wavelengths such that the total power is
% model.MC.P (which has a default value of 1 W). 
%
% When simulating broadband sources, you must specify model.MC.wavelength
% (and/or model.FMC.wavelength) as a 1D array of wavelengths. The specified
% simulationTimeRequested or nPhotonsRequested will be split between all
% the wavelengths.
% 
% For fluorescence light, the individual media's output spectra are
% determined by the ES (emission spectrum) specified in the mediaProperties
% function, although you must still choose which wavelengths to simulate in
% the model.FMC.wavelength array. For fluorescence emission, the total
% power of emission for each fluorescing voxel is determined by its
% medium's quantum yield. The emission spectrum only determines how that
% output power is distributed among the wavelengths.

%% MCmatlab abbreviations
% G: Geometry, MC: Monte Carlo, FMC: Fluorescence Monte Carlo, HS: Heat
% simulation, M: Media array, FR: Fluence rate, FD: Fractional damage.
%
% There are also some optional abbreviations you can use when referencing
% object/variable names: LS = lightSource, LC = lightCollector, FPID =
% focalPlaneIntensityDistribution, AID = angularIntensityDistribution, NI =
% normalizedIrradiance, NFR = normalizedFluenceRate.
%
% For example, "model.MC.LS.FPID.radialDistr" is the same as
% "model.MC.lightSource.focalPlaneIntensityDistribution.radialDistr"

%% Geometry definition
MCmatlab.closeMCmatlabFigures();
model = MCmatlab.model;

model.G.nx                = 100; % Number of bins in the x direction
model.G.ny                = 100; % Number of bins in the y direction
model.G.nz                = 100; % Number of bins in the z direction
model.G.Lx                = .1; % [cm] x size of simulation cuboid
model.G.Ly                = .1; % [cm] y size of simulation cuboid
model.G.Lz                = .1; % [cm] z size of simulation cuboid

model.G.mediaPropertiesFunc = @mediaPropertiesFunc; % Media properties defined as a function at the end of this file
model.G.geomFunc          = @geometryDefinition; % Function to use for defining the distribution of media in the cuboid. Defined at the end of this m file.

model = plot(model,'G');

%% Monte Carlo simulation
model.MC.simulationTimeRequested  = .1; % [min] Time duration of the simulation

model.MC.matchedInterfaces        = true; % Assumes all refractive indices are the same
model.MC.boundaryType             = 1; % 0: No escaping boundaries, 1: All cuboid boundaries are escaping, 2: Top cuboid boundary only is escaping, 3: Top and bottom boundaries are escaping, while the side boundaries are cyclic
model.MC.wavelength               = linspace(300,500,11); % [nm] Array of wavelength(s) for which to run incident-light Monte Carlo simulations
model.MC.spectrumFunc             = @Sfunc; % Defined just above the geometry function, later in this file

model.MC.lightSource.sourceType   = 4; % 0: Pencil beam, 1: Isotropically emitting line or point source, 2: Infinite plane wave, 3: Laguerre-Gaussian LG01 beam, 4: Radial-factorizable beam (e.g., a Gaussian beam), 5: X/Y factorizable beam (e.g., a rectangular LED emitter)
model.MC.lightSource.focalPlaneIntensityDistribution.radialDistr = 1; % Radial focal plane intensity distribution - 0: Top-hat, 1: Gaussian, Array: Custom. Doesn't need to be normalized.
model.MC.lightSource.focalPlaneIntensityDistribution.radialWidth = .03; % [cm] Radial focal plane 1/e^2 radius if top-hat or Gaussian or half-width of the full distribution if custom
model.MC.lightSource.angularIntensityDistribution.radialDistr = 0; % Radial angular intensity distribution - 0: Top-hat, 1: Gaussian, 2: Cosine (Lambertian), Array: Custom. Doesn't need to be normalized.
model.MC.lightSource.angularIntensityDistribution.radialWidth = 0; % [rad] Radial angular 1/e^2 half-angle if top-hat or Gaussian or half-angle of the full distribution if custom. For a diffraction limited Gaussian beam, this should be set to model.MC.wavelength*1e-9/(pi*model.MC.lightSource.focalPlaneIntensityDistribution.radialWidth*1e-2))
model.MC.lightSource.xFocus       = 0; % [cm] x position of focus
model.MC.lightSource.yFocus       = 0; % [cm] y position of focus
model.MC.lightSource.zFocus       = 0; % [cm] z position of focus
model.MC.lightSource.theta        = 0; % [rad] Polar angle of beam center axis
model.MC.lightSource.phi          = 0; % [rad] Azimuthal angle of beam center axis


model = runMonteCarlo(model);
model = plot(model,'MC');

%% Fluorescence Monte Carlo simulation
model.FMC.simulationTimeRequested  = .1; % [min] Time duration of the simulation

model.FMC.matchedInterfaces        = true; % Assumes all refractive indices are the same
model.FMC.boundaryType             = 1; % 0: No escaping boundaries, 1: All cuboid boundaries are escaping, 2: Top cuboid boundary only is escaping, 3: Top and bottom boundaries are escaping, while the side boundaries are cyclic
model.FMC.wavelength               = linspace(500,800,11); % [nm] Array of wavelength(s) for which to run fluorescence-light Monte Carlo simulations


model = runMonteCarlo(model,'fluorescence');
model = plot(model,'FMC');

%% Spectrum function
% The spectral power for the excitation light:
function spectralpower = Sfunc(wavelength)
  spectralpower = exp(-(wavelength-400).^2/(50).^2); % A simple Gaussian
end

%% Geometry function(s) (see readme for details)
function M = geometryDefinition(X,Y,Z,parameters)
  % Blood vessel example:
  zsurf = 0.02;
  epd_thick = 0.006;
  vesselradius  = 0.0100;
  vesseldepth = 0.04;
  M = ones(size(X)); % fill background with air
  M(Z > zsurf) = 2; % epidermis
  M(Z > zsurf + epd_thick) = 3; % dermis
  M(X.^2 + (Z - (zsurf + vesseldepth)).^2 < vesselradius^2) = 4; % blood
  R_sphere = 0.01;
  M(X.^2 + (Y - 0.03).^2 + (Z - 0.038).^2 < R_sphere^2) = 5; % stained cancerous tissue
  M(X.^2 + (Y + 0.03).^2 + (Z - zsurf).^2 < R_sphere^2 & Z < zsurf) = 6; % Lucifer yellow fluorescer
end

%% Media Properties function (see readme for details)
function mediaProperties = mediaPropertiesFunc(parameters)
  mediaProperties = MCmatlab.mediumProperties;

  j=1;
  mediaProperties(j).name  = 'air';
  mediaProperties(j).mua   = 1e-8; % [cm^-1]
  mediaProperties(j).mus   = 1e-8; % [cm^-1]
  mediaProperties(j).g     = 1.0;

  j=2;
  mediaProperties(j).name  = 'epidermis';
  mediaProperties(j).mua = @func_mua2;
  function mua = func_mua2(wavelength)
    B = 0; % Blood content
    S = 0.75; % Blood oxygen saturation
    W = 0.75; % Water content
    M = 0.03; % Melanin content
    F = 0; % Fat content
    mua = calc_mua(wavelength,S,B,W,F,M); % Jacques "Optical properties of biological tissues: a review" eq. 12
  end

  mediaProperties(j).mus = @func_mus2;
  function mus = func_mus2(wavelength)
    aPrime = 40; % musPrime at 500 nm
    fRay = 0; % Fraction of scattering due to Rayleigh scattering
    bMie = 1; % Scattering power for Mie scattering
    g = 0.9; % Scattering anisotropy
    mus = calc_mus(wavelength,aPrime,fRay,bMie,g); % Jacques "Optical properties of biological tissues: a review" eq. 2
  end
  mediaProperties(j).g   = 0.9;

  j=3;
  mediaProperties(j).name = 'dermis';
  mediaProperties(j).mua = @func_mua3;
  function mua = func_mua3(wavelength)
    B = 0.002; % Blood content
    S = 0.67; % Blood oxygen saturation
    W = 0.65; % Water content
    M = 0; % Melanin content
    F = 0; % Fat content
    mua = calc_mua(wavelength,S,B,W,F,M); % Jacques "Optical properties of biological tissues: a review" eq. 12
  end

  mediaProperties(j).mus = @func_mus3;
  function mus = func_mus3(wavelength)
    aPrime = 42.4; % musPrime at 500 nm
    fRay = 0.62; % Fraction of scattering due to Rayleigh scattering
    bMie = 1; % Scattering power for Mie scattering
    g = 0.9; % Scattering anisotropy
    mus = calc_mus(wavelength,aPrime,fRay,bMie,g); % Jacques "Optical properties of biological tissues: a review" eq. 2
  end
  mediaProperties(j).g   = 0.9;

  j=4;
  mediaProperties(j).name  = 'blood';
  mediaProperties(j).mua = @func_mua4;
  function mua = func_mua4(wavelength)
    B = 1; % Blood content
    S = 0.75; % Blood oxygen saturation
    W = 0.95; % Water content
    M = 0; % Melanin content
    F = 0; % Fat content
    mua = calc_mua(wavelength,S,B,W,F,M); % Jacques "Optical properties of biological tissues: a review" eq. 12
  end
  
  mediaProperties(j).mus = @func_mus4;
  function mus = func_mus4(wavelength)
    aPrime = 10; % musPrime at 500 nm
    fRay = 0; % Fraction of scattering due to Rayleigh scattering
    bMie = 1; % Scattering power for Mie scattering
    g = 0.9; % Scattering anisotropy
    mus = calc_mus(wavelength,aPrime,fRay,bMie,g); % Jacques "Optical properties of biological tissues: a review" eq. 2
  end
  mediaProperties(j).g   = 0.9;

  j=5;
  mediaProperties(j).name  = 'stained cancerous tissue';
  mediaProperties(j).mua = @func_mua5;
  function mua = func_mua5(wavelength)
    % We model the stained cancerous tissue absorption as a sum of the
    % absorption of dermis and an extra contribution from absorption around
    % 450 nm due to the fluorophore:
    B = 0.002; % Blood content
    S = 0.67; % Blood oxygen saturation
    W = 0.65; % Water content
    M = 0; % Melanin content
    F = 0; % Fat content
    mua = calc_mua(wavelength,S,B,W,F,M) + ... % Jacques "Optical properties of biological tissues: a review" eq. 12
          100*exp(-(wavelength-450).^2/(50).^2); % A simple Gaussian
  end

  mediaProperties(j).mus = @func_mus5;
  function mus = func_mus5(wavelength)
    aPrime = 42.4; % musPrime at 500 nm
    fRay = 0.62; % Fraction of scattering due to Rayleigh scattering
    bMie = 1; % Scattering power for Mie scattering
    g = 0.9; % Scattering anisotropy
    mus = calc_mus(wavelength,aPrime,fRay,bMie,g); % Jacques "Optical properties of biological tissues: a review" eq. 2
  end
  mediaProperties(j).g   = 0.9;

  mediaProperties(j).QY   = @QYfunc5; % Fluorescence quantum yield
  function QY = QYfunc5(wavelength)
    QY = 0.1*exp(-(wavelength-450).^2/(50).^2); % A simple Gaussian
  end
  mediaProperties(j).ES = @ESfunc5;
  function ES = ESfunc5(wavelength)
    ES = exp(-(wavelength-650).^2/(50).^2); % A simple Gaussian
  end

  j=6;
  mediaProperties(j).name  = 'Lucifer Yellow CH in water';
  mediaProperties(j).mua = @muafunc6;
  function mua = muafunc6(wavelength)
    absorption = readtable("helperfuncs/LuciferYellowCHinWater-abs.txt"); % The absorption spectrum was downloaded from omlc.org
    mua = interp1(absorption.Wavelength_nm,absorption.MolarExtinction_percmperM,wavelength,'linear',1e-8)/100; % interpolate from the raw data to the wavelength in question. Assuming a concentration of 10 mM. Allow extrapolation with value 1e-8 (mua may not be zero).
  end
  mediaProperties(j).mus = 10; % [cm^-1] This number is arbitrarily chosen in this case.
  mediaProperties(j).g   = 0.9;
  mediaProperties(j).QY   = 0.21; % Fluorescence quantum yield
  mediaProperties(j).ES = @ESfunc6;
  function ES = ESfunc6(wavelength)
    emission = readtable("helperfuncs/LuciferYellowCHinWater-ems.txt"); % The emission spectrum was downloaded from omlc.org
    ES = interp1(emission.Wavelength_nm,emission.Emission_AU,wavelength,'linear',0); % interpolate from the raw data to the wavelength in question. This doesn't need to be scaled, as MCmatlab normalises this function internally.
  end

end