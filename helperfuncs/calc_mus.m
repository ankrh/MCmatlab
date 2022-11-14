function mus = calc_mus(wavelength,aPrime,fRay,bMie,g)
fMie = 1 - fRay;
musPrime= aPrime*(fRay*(wavelength/500)^(-4) + fMie*(wavelength/500)^(-bMie)); % Jacques "Optical properties of biological tissues: a review" eq. 2
mus = musPrime/(1-g);
end