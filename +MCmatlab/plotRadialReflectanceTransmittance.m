function plotRadialReflectanceTransmittance(model)
  [X,Y] = ndgrid(model.G.x,model.G.y);
  R = sqrt(X.^2 + Y.^2);
  Runq = unique(R(:));
  Refl = NaN(size(Runq));
  Tran = NaN(size(Runq));
  for i = 1:numel(Runq)
    Refl(i) = mean(model.MC.NI_zneg(R == Runq(i)));
    Tran(i) = mean(model.MC.NI_zpos(R == Runq(i)));
  end
  Refl = [flipud(Refl(Runq ~= 0)) ; Refl];
  Tran = [flipud(Tran(Runq ~= 0)) ; Tran];
  Runq = [-flipud(Runq(Runq ~= 0)) ; Runq];
  MCmatlab.NdimSliderPlot(Refl,'nFig',43,'axisLabels',{'r [cm]','Reflectance [cm^{-2}]'},'axisValues',{Runq});
  MCmatlab.NdimSliderPlot(Tran,'nFig',44,'axisLabels',{'r [cm]','Transmittance [cm^{-2}]'},'axisValues',{Runq});
end