%% Decription
% This example uses the media and geometry of example 4 and serves to
% illustrate how the deposition criteria minMediumIdxToConsider and
% maxMediumIdxToConsider work. See also examples 8 and 19.
% 
% minMediumIdxToConsider and maxMediumIdxToConsider specify the minimum and
% maximum medium index that the rest of the deposition criteria apply to.
% For example, if you've specified minScatterings = 2,
% minMediumIdxToConsider = 4 and maxMediumIdxToConsider = 5, then only
% photons that have experienced at least 2 scattering events in media 4
% and/or 5 will deposit their weight in the output arrays. It doesn't
% matter how many scattering events have happened in other media.
% 
% Interface transitions and refractions will be counted if the medium
% transitioned *into* has index between minMediumIdxToConsider and
% maxMediumIdxToConsider.
% 
% In this example, we set minInterfaceTransitions = 1,
% minMediumIdxToConsider = 4 and maxMediumIdxToConsider = 4 so that we only
% see deposition and example paths from the photons that have entered the
% blood medium (index 4) at least once. We don't want to see deposition
% happen before the photon enters the blood vessel, so we set
% model.MC.depositionCriteria.evaluateOnlyAtEndOfLife = false.
% 
% The effect of the deposition criteria is most clearly visible in figure
% 5, where it is clear that all the photon paths originate on the surface
% of the blood medium.
% 
% Keep in mind that you can rearrange the media in the media definition as
% required so that the set of media that you want to consider in deposition
% criteria are in a contiguous interval.

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
model.MC.nExamplePaths = 100;

model.MC.matchedInterfaces        = true; % Assumes all refractive indices are the same
model.MC.boundaryType             = 1; % 0: No escaping boundaries, 1: All cuboid boundaries are escaping, 2: Top cuboid boundary only is escaping, 3: Top and bottom boundaries are escaping, while the side boundaries are cyclic
model.MC.wavelength               = 532; % [nm] Excitation wavelength, used for determination of optical properties for excitation light

model.MC.lightSource.sourceType   = 4; % 0: Pencil beam, 1: Isotropically emitting line or point source, 2: Infinite plane wave, 3: Laguerre-Gaussian LG01 beam, 4: Radial-factorizable beam (e.g., a Gaussian beam), 5: X/Y factorizable beam (e.g., a rectangular LED emitter)
model.MC.lightSource.focalPlaneIntensityDistribution.radialDistr = 0; % Radial focal plane intensity distribution - 0: Top-hat, 1: Gaussian, Array: Custom. Doesn't need to be normalized.
model.MC.lightSource.focalPlaneIntensityDistribution.radialWidth = .03; % [cm] Radial focal plane 1/e^2 radius if top-hat or Gaussian or half-width of the full distribution if custom
model.MC.lightSource.angularIntensityDistribution.radialDistr = 0; % Radial angular intensity distribution - 0: Top-hat, 1: Gaussian, 2: Cosine (Lambertian), Array: Custom. Doesn't need to be normalized.
model.MC.lightSource.angularIntensityDistribution.radialWidth = 0; % [rad] Radial angular 1/e^2 half-angle if top-hat or Gaussian or half-angle of the full distribution if custom. For a diffraction limited Gaussian beam, this should be set to model.MC.wavelength*1e-9/(pi*model.MC.lightSource.focalPlaneIntensityDistribution.radialWidth*1e-2))
model.MC.lightSource.xFocus       = 0; % [cm] x position of focus
model.MC.lightSource.yFocus       = 0; % [cm] y position of focus
model.MC.lightSource.zFocus       = 0; % [cm] z position of focus
model.MC.lightSource.theta        = 0; % [rad] Polar angle of beam center axis
model.MC.lightSource.phi          = 0; % [rad] Azimuthal angle of beam center axis

model.MC.depositionCriteria.minInterfaceTransitions = 1;
model.MC.depositionCriteria.minMediumIdxToConsider = 4;
model.MC.depositionCriteria.maxMediumIdxToConsider = 4;
model.MC.depositionCriteria.evaluateOnlyAtEndOfLife = false;

model = runMonteCarlo(model);
model = plot(model,'MC');

%% Geometry function(s) (see readme for details)
function M = geometryDefinition(X,Y,Z,parameters)
  % Blood vessel example:
  zsurf = 0.01;
  epd_thick = 0.006;
  vesselradius  = 0.0100;
  vesseldepth = 0.04;
  M = ones(size(X)); % fill background with water (gel)
  M(Z > zsurf) = 2; % epidermis
  M(Z > zsurf + epd_thick) = 3; % dermis
  M(X.^2 + (Z - (zsurf + vesseldepth)).^2 < vesselradius^2) = 4; % blood
end

%% Media Properties function (see readme for details)
function mediaProperties = mediaPropertiesFunc(parameters)
mediaProperties = MCmatlab.mediumProperties;

  j=1;
  mediaProperties(j).name  = 'water';
  mediaProperties(j).mua   = 0.00036; % [cm^-1]
  mediaProperties(j).mus   = 10; % [cm^-1]
  mediaProperties(j).g     = 1.0;
  mediaProperties(j).VHC   = 4.19; % [J cm^-3 K^-1]
  mediaProperties(j).TC    = 5.8e-3; % [W cm^-1 K^-1]

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
  mediaProperties(j).VHC = 3391*1.109e-3; % [J cm^-3 K^-1]
  mediaProperties(j).TC  = 0.37e-2; % [W cm^-1 K^-1]

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
  mediaProperties(j).VHC = 3391*1.109e-3; % [J cm^-3 K^-1]
  mediaProperties(j).TC  = 0.37e-2; % [W cm^-1 K^-1]

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
  mediaProperties(j).VHC = 3617*1.050e-3; % [J cm^-3 K^-1]
  mediaProperties(j).TC  = 0.52e-2; % [W cm^-1 K^-1]
  mediaProperties(j).E   = 422.5e3; % J/mol    PLACEHOLDER DATA ONLY
  mediaProperties(j).A   = 7.6e66; % 1/s        PLACEHOLDER DATA ONLY
end
