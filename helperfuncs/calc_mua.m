function mua = calc_mua(wavelength,S,B,W,F,M)
  persistent muadeoxy muafat muamel muaoxy muawater musp nmLIB
  if isempty(muadeoxy)
    load('spectralLIB.mat','muadeoxy','muafat','muamel','muaoxy','muawater','musp','nmLIB'); %#ok<NASGU>
  end
  mua_deoxy = interp1(nmLIB,muadeoxy,wavelength);
  mua_fat = interp1(nmLIB,muafat,wavelength);
  mua_mel = interp1(nmLIB,muamel,wavelength);
  mua_oxy = interp1(nmLIB,muaoxy,wavelength);
  mua_water = interp1(nmLIB,muawater,wavelength);

  % Jacques "Optical properties of biological tissues: a review" eq. 12:
  mua = B*S*mua_oxy + B*(1-S)*mua_deoxy + W*mua_water + F*mua_fat + M*mua_mel; % Equation 12 without the bilirubin and beta-Carotene terms
end